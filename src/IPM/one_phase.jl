function one_phase_solve(m::JuMP.Model)
    nlp_raw = MathProgNLPModel(m);
    return one_phase_solve(nlp_raw)
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
    intial_it = init(nlp, pars, timer);
    pause_advanced_timer(timer, "INIT")

    pause_advanced_timer(timer)
    print_timer_stats(timer)

    start_advanced_timer(timer)

    @assert(is_feasible(intial_it, pars.ls.comp_feas))
    iter, status, hist, t, err = one_phase_IPM(intial_it, pars, timer);

    pause_advanced_timer(timer)

    print_timer_stats(timer)

    return iter, status, hist, t, err
end

function switching_condition(iter::Class_iterate, last_step_was_superlinear::Bool, pars::Class_parameters)
    is_feas = is_feasible(iter, pars.ls.comp_feas_agg)
    dual_avg = scaled_dual_feas(iter, pars)

    if pars.primal_bounds_dual_feas
      dual_progress = dual_avg < norm(get_primal_res(iter), Inf)
    else
      dual_progress = dual_avg < get_mu(iter)
    end
    delta_small = get_delta(iter) < sqrt(get_mu(iter)) * max(0.1, norm(get_y(iter),Inf))
    lag_grad = norm(eval_grad_lag(iter,get_mu(iter)),1) < sum(iter.point.s .* iter.point.y) + norm(get_grad(iter) + iter.point.mu * eval_grad_r(iter),1) # + norm(get_primal_res(iter), Inf) + 1.0 #+ sqrt(norm(get_y(iter),Inf))

    be_aggressive = is_feas && dual_progress && lag_grad
    be_aggressive |= last_step_was_superlinear && dual_progress && lag_grad

    return be_aggressive
end

function one_phase_IPM(iter::Class_iterate, pars::Class_parameters, timer::class_advanced_timer)
  t = 0;
  progress = Array{alg_history,1}();
  filter = Array{Class_filter,1}();

      update!(iter, timer, pars) # is this necessary ????

      kkt_solver = pick_KKT_solver(pars);
      initialize!(kkt_solver, iter)

      head_progress()

      record_progress_first_it!(iter, kkt_solver, progress, pars)

      init_step_size = 1.0
      status = :success
      agg_next_step = false;
      dir_size_agg = Inf
      dir_size_stable = 0.0
      ls_info = false

      start_time = time()
      last_step_was_superlinear = false

      for t = 1:pars.term.max_it
             @assert(is_feasible(iter, pars.ls.comp_feas))

             for i = 1:pars.max_it_corrections
               tot_num_fac = 0; inertia_num_fac = 0;

               start_advanced_timer(timer, "misc/terminate")
               status = terminate(iter, pars)
               pause_advanced_timer(timer, "misc/terminate")

               if status != false
                 println("Terminated with ", status)
                 return iter, status, progress, t, false
               end

               if time() - start_time > pars.term.max_time
                 println("Terminated due to timeout")
                 return iter, :MAX_TIME, progress, t, false
               end

               start_advanced_timer(timer, "misc/checks")
               be_aggressive = switching_condition(iter, last_step_was_superlinear, pars)
               last_step_was_superlinear = false
               pause_advanced_timer(timer, "misc/checks")

               if i == 1
                   update_H!(iter, timer, pars)
                   @assert(is_updated(iter.cache))

                   form_system!(kkt_solver, iter, timer)

                   start_advanced_timer(timer, "ipopt_strategy")
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
                     println(pd("**"), pd("status"), pd("delta"), pd("step"))
                   end

                   for k = 1:100
                     #status, new_iter, ls_info = take_step!(iter, reduct_factors, kkt_solver, ls_mode, filter, pars, actual_min_step_size, timer)
                     status, new_iter, ls_info, reduct_factors = take_step2!(be_aggressive, iter, kkt_solver, filter, pars, timer)


                     if pars.output_level >= 6
                       println(pd("**"), pd(status), rd(get_delta(iter)), rd(ls_info.step_size_P), rd(norm(kkt_solver.dir.x,Inf)), rd(norm(kkt_solver.dir.y,Inf)), rd(norm(kkt_solver.dir.s,Inf)))
                     end

                     if status == :success
                       break
                     elseif i < 100 && get_delta(iter) < pars.delta.max
                       if pars.test.response_to_failure == :lag_delta_inc
                         set_delta(iter, max(norm(eval_grad_lag(iter,iter.point.mu),Inf) / norm(kkt_solver.dir.x,Inf),get_delta(iter) * pars.delta.inc, max(pars.delta.start, old_delta * pars.delta.dec)))
                       elseif pars.test.response_to_failure == :default
                         set_delta(iter, max(get_delta(iter) * pars.delta.inc, max(pars.delta.start, old_delta * pars.delta.dec)))
                       else
                          error("pars.test.response_to_failure parameter incorrectly set")
                       end
                       inertia = factor!(kkt_solver, get_delta(iter), timer)
                       tot_num_fac += 1
                     else
                       pause_advanced_timer(timer, "STEP/first")
                       pause_advanced_timer(timer, "STEP")
                       println("Terminated due to max delta while attempting to take step")
                       println("delta=$(get_delta(iter)), be_aggressive=$be_aggressive, status=$status")
                       println("dx = $(norm(kkt_solver.dir.x,2)), dy = $(norm(kkt_solver.dir.y,2)), ds = $(norm(kkt_solver.dir.s,2))")
                       @show reduct_factors #, ls_mode
                       @show ls_info
                       return iter, :MAX_DELTA, progress, t, false
                     end
                   end

                   pause_advanced_timer(timer, "STEP/first")
                   pause_advanced_timer(timer, "STEP")
               else
                   start_advanced_timer(timer, "STEP")
                   start_advanced_timer(timer, "STEP/correction")

                   status, new_iter, ls_info, reduct_factors = take_step2!(be_aggressive, iter, kkt_solver, filter, pars, timer)

                   if pars.superlinear_theory_mode && be_aggressive
                     if get_mu(new_iter) < get_mu(iter) * 0.1
                        last_step_was_superlinear = true
                     end
                   end

                   pause_advanced_timer(timer, "STEP/correction")
                   pause_advanced_timer(timer, "STEP")
               end

               if status == :success
                 iter = new_iter
                 if be_aggressive
                   dir_size_agg = norm(kkt_solver.dir.x, 2)
                 end
               end

               add!(filter, iter, pars)
               start_advanced_timer(timer, "misc/record_progress")
               record_progress!(t, be_aggressive ? "agg" : "stb", iter, kkt_solver, ls_info, reduct_factors, progress, inertia_num_fac, tot_num_fac, pars)
               pause_advanced_timer(timer, "misc/record_progress")
               @assert(is_updated_correction(iter.cache))
               check_for_nan(iter.point)

               if status != :success
                 break
               end
             end
      end

  println("Terminated due to max iterations reached")
  return iter, :MAX_IT, progress, pars.term.max_it, false
end
