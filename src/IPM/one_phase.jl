function one_phase_solve(m::JuMP.Model)
    nlp_raw = MathOptNLPModel(m);
    return one_phase_solve(nlp_raw)
end

function one_phase_solve(m::JuMP.Model, pars::Class_parameters)
    nlp_raw = MathOptNLPModel(m);
    return one_phase_solve(nlp_raw, pars::Class_parameters)
end


function one_phase_solve(nlp_raw::NLPModels.AbstractNLPModel)
    pars = Class_parameters()
    return one_phase_solve(nlp_raw, pars)
end

function one_phase_solve(nlp_raw::NLPModels.AbstractNLPModel, pars::Class_parameters)
    nlp = Class_CUTEst(nlp_raw);
    timer = class_advanced_timer()
    start_advanced_timer(timer)

    start_advanced_timer(timer, "INIT")
    if pars.init.init_style == :gertz
        intial_it = gertz_init(nlp, pars, timer); # Gertz, Michael, Jorge Nocedal, and A. Sartenar. "A starting point strategy for nonlinear interior methods." Applied mathematics letters 17.8 (2004): 945-952.
    elseif pars.init.init_style == :mehrotra
        intial_it = mehrotra_init(nlp, pars, timer);
    elseif pars.init.init_style == :LP
        intial_it = LP_init(nlp, pars, timer);
    else
        error("Init strategy does not exist")
    end
    pause_advanced_timer(timer, "INIT")

    pause_advanced_timer(timer)
    if pars.output_level >= 4
        print_timer_stats(timer)
    end

    start_advanced_timer(timer)

    @assert(is_feasible(intial_it, pars.ls.comp_feas))
    iter, status, hist, t, err = one_phase_IPM(intial_it, pars, timer);

    pause_advanced_timer(timer)

    if pars.output_level >= 3
        print_timer_stats(timer)
    end
    return iter, status, hist, t, err, timer
end

function switching_condition(iter::Class_iterate, last_step_was_superlinear::Bool, pars::Class_parameters)
    # should we take an aggressive step or not?
    is_feas = is_feasible(iter, pars.ls.comp_feas_agg)
    dual_avg = scaled_dual_feas(iter, pars)

    if pars.primal_bounds_dual_feas
      dual_progress = dual_avg < pars.aggressive_dual_threshold * LinearAlgebra.norm(get_primal_res(iter), Inf)
    else
      dual_progress = dual_avg < pars.aggressive_dual_threshold * get_mu(iter)
    end
    delta_small = get_delta(iter) < sqrt(get_mu(iter)) * max(0.1, LinearAlgebra.norm(get_y(iter),Inf))
    lag_grad = LinearAlgebra.norm(eval_grad_lag(iter,get_mu(iter)),1) < sum(iter.point.s .* iter.point.y) + LinearAlgebra.norm(get_grad(iter) + iter.point.mu * eval_grad_r(iter),1) # + LinearAlgebra.norm(get_primal_res(iter), Inf) + 1.0 #+ sqrt(LinearAlgebra.norm(get_y(iter),Inf))

    be_aggressive = is_feas && dual_progress && lag_grad
    be_aggressive |= last_step_was_superlinear && dual_progress && lag_grad

    return be_aggressive
end

function one_phase_IPM(iter::Class_iterate, pars::Class_parameters, timer::class_advanced_timer)
  #####################################################################
  # THE MAIN ALGORITHM
  # input:
  # iter = starting point
  # pars = parameters for running the algorithm
  # timer = code to time the algorithm
  #####################################################################

      # intialize code
      t = 0;
      rpt = 0.0
      progress = Array{alg_history2,1}();
      filter = Array{Class_filter,1}();

      update!(iter, timer, pars) # is this necessary ????

      kkt_solver = pick_KKT_solver(pars);
      initialize!(kkt_solver, iter)

      if pars.output_level >= 1
          head_progress()
      end

      display = pars.output_level >= 1
      record_progress_first_it!(progress, iter, kkt_solver, pars, display)
      if pars.output_level >= 4
          println("")
      end

      init_step_size = 1.0
      status = :success
      agg_next_step = false;
      dir_size_agg = Inf
      dir_size_stable = 0.0
      ls_info = false

      start_time = time()
      last_step_was_superlinear = false
      scale_update = false

      # check termination critieron at starting point
      start_advanced_timer(timer, "misc/terminate")
      status = terminate(iter, pars)
      pause_advanced_timer(timer, "misc/terminate")

      if status != false
        if pars.output_level >= 1
            println("Terminated with ", status)
        end
        return iter, status, progress, t, false
      end

      if time() - start_time > pars.term.max_time
          if pars.output_level >= 1
            println("Terminated due to timeout")
        end
        return iter, :MAX_TIME, progress, t, false
      end

      # run the algorithm
      for t = 1:pars.term.max_it
             @assert(is_feasible(iter, pars.ls.comp_feas))

             for i = 1:pars.max_it_corrections
                         if pars.output_level >= 4
                             println("================================== ITERATION $t, MINOR ITERATION $i ======================================")
                         end

                         if pars.output_level >= 5
                             println("Strict comp = ", maximum(min.(iter.point.s/LinearAlgebra.norm(iter.point.s,Inf),iter.point.y/LinearAlgebra.norm(iter.point.y,Inf))))
                         end

                       tot_num_fac = 0; inertia_num_fac = 0;

                       start_advanced_timer(timer, "misc/checks")
                       be_aggressive = switching_condition(iter, last_step_was_superlinear, pars) # should we take an aggressive step or not?
                       last_step_was_superlinear = false

                       pause_advanced_timer(timer, "misc/checks")

                       if i == 1
                           ##########################################################################################
                           ##### first step of iteration, we need to compute lag hessian and do a factorization #####
                           ##########################################################################################
                           update_H!(iter, timer, pars)
                           @assert(is_updated(iter.cache))

                           # form matrix that we are going to factorize
                           form_system!(kkt_solver, iter, timer)

                           start_advanced_timer(timer, "ipopt_strategy")
                           # Figures out what delta will give
                           fact_succeed, inertia_num_fac, new_delta = ipopt_strategy!(iter, kkt_solver, pars, timer)
                           tot_num_fac = inertia_num_fac
                           old_delta = get_delta(iter)
                           set_delta(iter, new_delta)
                           pause_advanced_timer(timer, "ipopt_strategy")

                           if fact_succeed != :success
                              return iter, :MAX_DELTA, progress, t, false
                           end

                           start_advanced_timer(timer, "STEP")
                           start_advanced_timer(timer, "STEP/first")

                           if pars.output_level >= 5
                             println(pd("**"), pd("status"), pd("delta"), pd("α_P"), pd("dx"), pd("dy"), pd("ds"))
                           end

                           for k = 1:100
                             #status, new_iter, ls_info = take_step!(iter, reduct_factors, kkt_solver, ls_mode, filter, pars, actual_min_step_size, timer)
                             step_status, new_iter, ls_info, reduct_factors = take_step2!(be_aggressive, iter, kkt_solver, filter, pars, timer)


                             if pars.output_level >= 6
                               println("")
                               println(pd("**"), pd(step_status), rd(get_delta(iter)), rd(ls_info.step_size_P), rd(LinearAlgebra.norm(kkt_solver.dir.x,Inf)), rd(LinearAlgebra.norm(kkt_solver.dir.y,Inf)), rd(LinearAlgebra.norm(kkt_solver.dir.s,Inf)))
                             end

                             if step_status == :success
                               break
                             elseif i < 100 && get_delta(iter) < pars.delta.max
                               if pars.test.response_to_failure == :lag_delta_inc
                                 set_delta(iter, max(LinearAlgebra.norm(eval_grad_lag(iter,iter.point.mu),Inf) / LinearAlgebra.norm(kkt_solver.dir.x,Inf),get_delta(iter) * pars.delta.inc, max(pars.delta.start, old_delta * pars.delta.dec)))
                               elseif pars.test.response_to_failure == :default
                                 set_delta(iter, max(get_delta(iter) * pars.delta.inc, max(pars.delta.start, old_delta * pars.delta.dec)))
                               else
                                  error("pars.test.response_to_failure parameter incorrectly set")
                               end
                               inertia = factor!(kkt_solver, get_delta(iter), timer)
                               tot_num_fac += 1
                             elseif LinearAlgebra.norm(comp(iter),Inf) > 1e-14
                                 warn("Error ... large delta causing issues")
                                iter.point.y = iter.point.mu ./ iter.point.s
                                step_status = :success
                                break
                             else
                               pause_advanced_timer(timer, "STEP/first")
                               pause_advanced_timer(timer, "STEP")
                               println("Terminated due to max delta while attempting to take step")
                               println("delta=$(get_delta(iter)), be_aggressive=$be_aggressive, status=$step_status")
                               println("dx = $(LinearAlgebra.norm(kkt_solver.dir.x,2)), dy = $(LinearAlgebra.norm(kkt_solver.dir.y,2)), ds = $(LinearAlgebra.norm(kkt_solver.dir.s,2))")
                               @show reduct_factors #, ls_mode
                               @show ls_info
                               return iter, :MAX_DELTA, progress, t, false
                             end
                           end

                           pause_advanced_timer(timer, "STEP/first")
                           pause_advanced_timer(timer, "STEP")
                       else
                           #######################################
                           ### corrections, reuse factorization ##
                           #######################################
                           start_advanced_timer(timer, "STEP")
                           start_advanced_timer(timer, "STEP/correction")

                           step_status, new_iter, ls_info, reduct_factors = take_step2!(be_aggressive, iter, kkt_solver, filter, pars, timer)

                           if pars.superlinear_theory_mode && be_aggressive
                             if get_mu(new_iter) < get_mu(iter) * 0.1
                                last_step_was_superlinear = true
                             end
                           end

                           pause_advanced_timer(timer, "STEP/correction")
                           pause_advanced_timer(timer, "STEP")
                       end

                       if step_status == :success
                         iter = new_iter
                         if be_aggressive
                           dir_size_agg = LinearAlgebra.norm(kkt_solver.dir.x, 2)
                         end
                       end

                       add!(filter, iter, pars) # update filter

                       start_advanced_timer(timer, "misc/terminate")
                       status = terminate(iter, pars) # check termination criteron
                       pause_advanced_timer(timer, "misc/terminate")

                       start_advanced_timer(timer, "misc/record_progress")
                       # output to the console
                       output_level = pars.output_level
                       display = output_level >= 4 || (output_level >= 3 && i == 1) || (output_level == 2 && t % 10 == 1 && i == 1) || (status != false && output_level >= 1)

                       record_progress!(progress, t, be_aggressive ? "agg" : "stb", iter, kkt_solver, ls_info, reduct_factors, inertia_num_fac, tot_num_fac, pars, display)
                       if pars.output_level >= 4
                           println("")
                       end

                       pause_advanced_timer(timer, "misc/record_progress")
                       @assert(is_updated_correction(iter.cache))
                       check_for_nan(iter.point)

                       # if termination criteron is satisfied stop the algorithm
                       if status != false
                         if pars.output_level >= 1
                             println("Terminated with ", status)
                         end
                         return iter, status, progress, t, false
                       end

                       if time() - start_time > pars.term.max_time
                         if pars.output_level >= 1
                             println("Terminated due to timeout")
                         end
                         return iter, :MAX_TIME, progress, t, false
                       end

                       if step_status != :success
                         break
                       end
             end
      end

  if pars.output_level >= 1
      println("Terminated due to max iterations reached")
  end
  return iter, :MAX_IT, progress, pars.term.max_it, false
end
