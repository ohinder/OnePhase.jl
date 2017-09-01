function primal_update_of_dual!(iter, pars)
    D1 = norm(eval_grad_lag(iter, iter.point.y), 2)
    y_primal = iter.point.mu ./ iter.point.s
    D2 = norm(eval_grad_lag(iter, y_primal), 2)
    if D2 < D1
      iter.point.y = y_primal
    end
end

function mu_func(iter)
    return get_mu(iter) #norm(get_primal_res(iter),Inf) #/(1 + sqrt(get_mu(iter))) # norm(get_grad(iter),Inf)
end

function one_phase_IPM(iter::Class_iterate, pars::Class_parameters, timer::class_advanced_timer)
  t = 0;
  progress = Array{alg_history,1}();
  filter = Array{Class_filter,1}();

  #try
      #iter.point.x += rand(length(iter.point.x)) / 10.0
      if pars.use_prox
        convert_to_prox!(iter, pars, mu_func(iter));
      end
      update!(iter, timer, pars)

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
            step_type = nothing; new_iter = nothing;

             @assert(is_feasible(iter, pars.comp_feas))

             for i = 1:pars.max_it_corrections
               status = terminate(iter, pars)
               if status != false
                 println("Terminated with ", status)
                 return iter, status, progress, t, false
               end

               if time() - start_time > pars.MAX_TIME
                 println("Terminated due to timeout")
                 return iter, :MAX_TIME, progress, t, false
               end

               primal_inf = norm(iter.primal_residual_intial, Inf) * iter.point.primal_scale

               feas_obj = max(0.0,-mean(iter.point.y .* iter.cache.cons))
               if pars.adaptive_mu == :none
                 mu_est = iter.point.mu
               elseif pars.adaptive_mu == :test5
                 mu_est = mean( (iter.point.y + 0.1) .* - get_primal_res(iter)) / 2.0
               elseif pars.adaptive_mu == :test6
                 mu_est = mean(iter.point.y .* iter.point.s) * 0.5 + feas_obj * 0.5
               elseif pars.adaptive_mu == :test7
                 mu_avg_guarded = min(iter.point.mu, mean(iter.point.y .* iter.point.s))
                 mu_est = max(mu_avg_guarded, min(10.0 * iter.point.mu, feas_obj)) + norm(get_primal_res(iter),Inf) * 0.01
               elseif pars.adaptive_mu == :test8
                 mu_est = mean(iter.point.y .* iter.point.s) / 10.0 + feas_obj + norm(get_primal_res(iter),Inf) * 0.01
               elseif pars.adaptive_mu == :test9
                 mu_est = mean(iter.point.y .* iter.point.s) + mean( (iter.point.y + 0.1) .* - get_primal_res(iter))
               elseif pars.adaptive_mu == :test10
                 mu_est = dot(iter.point.y, iter.point.s - get_primal_res(iter)) / (2.0 * length(iter.point.y))
               elseif pars.adaptive_mu == :test11
                 mu_est = mean(iter.point.y .* iter.point.s)
               elseif pars.adaptive_mu == :test12
                 mu_est = mean(iter.point.y .* iter.point.s) + feas_obj
               elseif pars.adaptive_mu == :paper
                 mu_est_unguarded = mean(iter.point.y .* iter.point.s) + feas_obj
                 Θ = norm(get_primal_res(iter),Inf)
                 mu_est = min(Θ * 1e5, max(mu_est_unguarded, Θ * 1e-5))
               elseif pars.adaptive_mu == :paper2
                 mu_est_unguarded = 0.9 * mean(iter.point.y .* iter.point.s) + feas_obj * 0.2
                 Θ = norm(get_primal_res(iter),Inf)
                 mu_est = min(Θ * 1e5, max(mu_est_unguarded, Θ * 1e-5))
               end

               is_feas = is_feasible(iter, pars.comp_feas_agg) && norm(comp(iter),Inf) < pars.comp_feas_agg_inf
               #dual_progress = scaled_dual_feas(iter, pars) < threshold
               #mu_est = iter.point.mu
               #mu_est = mean(iter.point.y .* iter.point.s) * 0.5 + feas_obj * 0.5
               #mu_est = max(min(100.0 * norm(get_primal_res(iter),Inf), iter.point.mu, mean(iter.point.y .* iter.point.s) ), feas_obj) + norm(get_primal_res(iter),Inf) * 0.01
               #mu_est = mean(iter.point.y .* iter.point.s)
               #mu_est = (mean(iter.point.y .* iter.point.s) + mean(-iter.point.y .* get_primal_res(iter))) / 4.0
               mu_ub = (mean(iter.point.y) + 1.0) * mean(-get_primal_res(iter))
               #dual_avg = norm(eval_grad_lag(iter),1) / length(iter.point.mu)
               dual_avg = scaled_dual_feas(iter, pars)

               if pars.primal_bounds_dual_feas
                 dual_progress = dual_avg < norm(get_primal_res(iter), Inf)
               else
                 dual_progress = dual_avg < get_mu(iter)
               end
               # * 10.0
               #dual_progress = dual_avg < mu_ub * 10.0
               delta_small = get_delta(iter) < sqrt(get_mu(iter)) * (1.0 + norm(get_y(iter),Inf))
               #
               lag_grad = norm(eval_grad_lag(iter),Inf) < get_mu(iter) + norm(get_grad(iter),Inf) + (norm(get_primal_res(iter), Inf) + 1.0) #+ sqrt(norm(get_y(iter),Inf))
               #@show norm(get_grad(iter),Inf)

               no_stall = true #(ls_info == false || ls_info.step_size_P > 0.1)
               be_aggressive = no_stall && is_feas && (lag_grad || !pars.lag_grad_test) && dual_progress && (delta_small || !pars.inertia_test)


                 if be_aggressive
                     if norm(get_primal_res(iter),Inf) > pars.tol || !pars.pause_primal
                       step_type = "agg"
                       reduct_factors = pars.aggressive_reduct_factors #Reduct_affine()
                     else
                       reduct_factors = Class_reduction_factors(1.0, 0.0, 0.0)
                       step_type = "mu"
                     end
                     #reduct_factors = Class_reduction_factors(0.5, 0.5, 0.5)
                     #ls_mode = pars.ls_mode_agg;

                     q = pars.min_step_size_agg_ratio * min(1.0, 1.0 / maximum(- get_primal_res(iter) ./ iter.point.s))
                     actual_min_step_size = q
                     #@show q
                     iter.point.mu = mu_est
                     centre_dual!(iter.point, (pars.comp_feas_agg + pars.comp_feas)/2.0) # ??????
                 else
                      reduct_factors = pars.stable_reduct_factors #Reduct_stable()
                      step_type = "stb"
                      actual_min_step_size = pars.min_step_size_stable
                 end

               if i == 1
                   #primal_update_of_dual!(iter, pars)
                   update_H!(iter, timer, pars)
                   @assert(is_updated(iter.cache))

                   ipopt_strategy!(iter, kkt_solver, pars, timer)

                   start_advanced_timer(timer, "STEP/first")

                   if pars.output_level >= 5
                     println(pd("**"), pd("status"), pd("delta"), pd("step"))
                   end

                   for i = 1:100
                     if be_aggressive
                       ls_mode = pars.ls_mode_agg;
                       if pars.use_prox
                         update_prox!(iter, pars, 0.0)
                         update!(iter, timer, pars)
                       end
                     else
                       if pars.use_prox
                         update_prox!(iter, pars, mu_func(iter))
                         update!(iter, timer, pars)
                       end

                       if get_delta(iter) == 0.0
                         ls_mode = pars.ls_mode_stable_delta_zero
                         #pars.ls_mode_stable_correction;
                         #;
                       else
                         ls_mode = pars.ls_mode_stable_correction;
                       end
                     end

                     status, new_iter, ls_info = take_step!(iter, reduct_factors, kkt_solver, ls_mode, filter, pars, actual_min_step_size, timer)

                     if pars.output_level >= 6
                       println(pd("**"), pd(status), rd(get_delta(iter)), rd(ls_info.step_size_P), rd(norm(kkt_solver.dir.x,Inf)), rd(norm(kkt_solver.dir.y,Inf)), rd(norm(kkt_solver.dir.s,Inf)))
                     end

                     if status == :success
                       #@assert(norm(iter.point.x - new_iter.point.x) > 0)
                       break
                     elseif i < 100
                       #8.0
                       #println("increase delta")
                       set_delta(iter, max(get_delta(iter) * 20.0, 1e-8))
                       inertia = factor!(kkt_solver, get_delta(iter), timer)
                     else
                       pause_advanced_timer(timer, "STEP/first")
                       println("Terminated due to max delta")
                       println("inertia=$(inertia)")
                       println("delta=$(get_delta(iter)), step_type=$step_type, min_step_size=$actual_min_step_size, status=$status")
                       println("dx = $(norm(kkt_solver.dir.x,2)), dy = $(norm(kkt_solver.dir.y,2)), ds = $(norm(kkt_solver.dir.s,2))")
                       @show reduct_factors, ls_mode
                       @show ls_info
                       return iter, :MAX_DELTA, progress, t, false
                     end
                   end

                   pause_advanced_timer(timer, "STEP/first")
               else
                   start_advanced_timer(timer, "STEP/correction")
                   if be_aggressive
                     ls_mode = pars.ls_mode_agg;
                   elseif get_delta(iter) == 0.0
                     #ls_mode = pars.ls_mode_stable_correction;
                     ls_mode = pars.ls_mode_stable_delta_zero;
                   else
                     ls_mode = pars.ls_mode_stable_correction;
                   end

                   status, new_iter, ls_info = take_step!(iter, reduct_factors, kkt_solver, ls_mode, filter, pars, actual_min_step_size, timer)

                   pause_advanced_timer(timer, "STEP/correction")
               end

               if status == :success
                 iter = new_iter
                 if step_type == "agg"
                   dir_size_agg = norm(kkt_solver.dir.x, 2)
                    if pars.use_prox
                       update_prox!(iter, pars, mu_func(iter))
                       update!(iter, timer, pars)
                    end
                 end
               end

               add!(filter, iter, pars)
               record_progress!(t,  step_type, iter, kkt_solver, ls_info, reduct_factors, progress, pars)
               @assert(is_updated_correction(iter.cache))
               check_for_nan(iter.point)

               if status != :success
                 break
               end
             end
      end
    #=
    catch(e)
      if isa(e, Eval_NaN_error)
        println("Terminated with Eval_NaN_error")
        return iter, :NaN_ERR, progress, t, e
      else
        println("Terminated with unknown error")
        return iter, :UNKNOWN_ERR, progress, t, e
      end
    end=#

  println("Terminated due to max iterations reached")
  return iter, :MAX_IT, progress, pars.max_it, false
end
