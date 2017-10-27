function primal_update_of_dual!(iter, pars)
    D1 = norm(eval_grad_lag(iter, iter.point.mu, iter.point.y), 2)
    y_primal = iter.point.mu ./ iter.point.s
    D2 = norm(eval_grad_lag(iter, iter.point.mu, y_primal), 2)
    if D2 < D1
      iter.point.y = y_primal
    end
end

function mu_func(iter)
    return get_mu(iter) #norm(get_primal_res(iter),Inf) #/(1 + sqrt(get_mu(iter))) # norm(get_grad(iter),Inf)
end

function one_phase_solve(m::JuMP.Model)
    nlp_raw = MathProgNLPModel(m);
    return one_phase_solve(nlp_raw)
end

function one_phase_solve(nlp_raw)
    nlp = Class_CUTEst(nlp_raw);
    timer = class_advanced_timer()
    start_advanced_timer(timer)
    #include("include.jl")
    #intial_it = initial_point_satisfy_bounds(nlp, my_par)
    start_advanced_timer(timer, "INIT")
    my_par = Class_parameters()
    intial_it = init(nlp, my_par, timer);
    #intial_it.point.mu *= 10.0
    #intial_it.point.y *= 10.0
    pause_advanced_timer(timer, "INIT")

    pause_advanced_timer(timer)
    print_timer_stats(timer)

    start_advanced_timer(timer)

    @assert(is_feasible(intial_it, my_par.comp_feas))
    iter, status, hist, t, err = one_phase_IPM(intial_it, my_par, timer);

    pause_advanced_timer(timer)

    print_timer_stats(timer)

    return iter
end

function one_phase_IPM(iter::Class_iterate, pars::Class_parameters, timer::class_advanced_timer)
  t = 0;
  progress = Array{alg_history,1}();
  filter = Array{Class_filter,1}();

  #try
      #iter.point.x += rand(length(iter.point.x)) / 10.0
      #if pars.use_prox
      #  convert_to_prox!(iter, pars, mu_func(iter));
      #end
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

      for t = 1:pars.max_it
             @assert(is_feasible(iter, pars.comp_feas))

             for i = 1:pars.max_it_corrections
               tot_num_fac = 0; inertia_num_fac = 0;

               start_advanced_timer(timer, "misc/terminate")
               status = terminate(iter, pars)
               pause_advanced_timer(timer, "misc/terminate")

               if status != false
                 println("Terminated with ", status)
                 return iter, status, progress, t, false
               end

               if time() - start_time > pars.MAX_TIME
                 println("Terminated due to timeout")
                 return iter, :MAX_TIME, progress, t, false
               end

               start_advanced_timer(timer, "misc/checks")
               primal_inf = norm(iter.primal_residual_intial, Inf) * iter.point.primal_scale

               mu_est = iter.point.mu

               is_feas = is_feasible(iter, pars.comp_feas_agg) && norm(comp(iter),Inf) < pars.comp_feas_agg_inf
               #dual_avg = norm(eval_grad_lag(iter),1) / length(iter.point.mu)
               dual_avg = scaled_dual_feas(iter, pars)

               if pars.primal_bounds_dual_feas
                 dual_progress = dual_avg < norm(get_primal_res(iter), Inf)
               else
                 dual_progress = dual_avg < get_mu(iter) #/ 10.0
               end
               # * 10.0
               delta_small = get_delta(iter) < sqrt(get_mu(iter)) * max(0.1, norm(get_y(iter),Inf))
               #
               #lag_grad = norm(eval_grad_lag(iter,get_mu(iter)),Inf) < get_mu(iter) + norm(get_grad(iter),Inf) + (norm(get_primal_res(iter), Inf) + 1.0) #+ sqrt(norm(get_y(iter),Inf))
               lag_grad = norm(eval_grad_lag(iter,get_mu(iter)),Inf) < max(get_mu(iter)/ pars.comp_feas_agg, norm(get_grad(iter),Inf)) # + norm(get_primal_res(iter), Inf) + 1.0 #+ sqrt(norm(get_y(iter),Inf))
               #lag_grad = norm(eval_grad_lag(iter,get_mu(iter)),Inf) < max(get_mu(iter), norm(get_grad(iter),Inf)) # + norm(get_primal_res(iter), Inf) + 1.0 #+ sqrt(norm(get_y(iter),Inf))
               #@show norm(get_grad(iter),Inf)
               #&& lag_grad
               be_aggressive = is_feas && dual_progress && lag_grad && (delta_small || !pars.inertia_test)
               #@show is_feas, dual_progress, lag_grad
               pause_advanced_timer(timer, "misc/checks")

               if i == 1
                   update_H!(iter, timer, pars)
                   @assert(is_updated(iter.cache))

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
                       #@assert(norm(iter.point.x - new_iter.point.x) > 0)
                       break
                     elseif i < 100 && get_delta(iter) < pars.delta.max
                       #8.0
                       #println("increase delta")
                       set_delta(iter, max(get_delta(iter) * pars.delta.inc, max(pars.delta.start, old_delta * pars.delta.dec)))
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

                   pause_advanced_timer(timer, "STEP/correction")
                   pause_advanced_timer(timer, "STEP")
               end

               if status == :success
                 iter = new_iter
                 if be_aggressive
                   dir_size_agg = norm(kkt_solver.dir.x, 2)
                    #if pars.use_prox
                    #   update_prox!(iter, pars, mu_func(iter))
                    #   update!(iter, timer, pars)
                    #end
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
  return iter, :MAX_IT, progress, pars.max_it, false
end
