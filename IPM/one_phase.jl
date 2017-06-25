

function one_phase_IPM(iter::Class_iterate, pars::Class_parameters, timer::class_advanced_timer)
  t = 0;
  progress = Array{alg_history,1}();
  filter = Array{Class_filter,1}();

  #try
      convert_to_prox!(iter, pars);

      kkt_solver = pick_KKT_solver(pars);
      initialize!(kkt_solver, iter)

      head_progress()

      record_progress_first_it!(iter, kkt_solver, progress, pars)

      init_step_size = 1.0
      success = :success
      agg_next_step = false;
      dir_size_agg = Inf
      dir_size_stable = 0.0

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

               #=if pars.threshold_type == :mu
                 threshold = pars.aggressive_dual_threshold * get_mu(iter)
               elseif pars.threshold_type == :mu_primal
                  threshold = pars.aggressive_dual_threshold * min(get_mu(iter),norm(get_primal_res(iter),Inf))
               elseif pars.threshold_type == :primal
                   threshold = pars.aggressive_dual_threshold * norm(get_primal_res(iter),Inf)
               else
                 error("threshold_type not known")
               end=#

               is_feas = is_feasible(iter, pars.comp_feas_agg) && norm(comp(iter),Inf) < pars.comp_feas_agg_inf
               #dual_progress = scaled_dual_feas(iter, pars) < threshold

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
               end

               #mu_est = iter.point.mu
               #mu_est = mean(iter.point.y .* iter.point.s) * 0.5 + feas_obj * 0.5
               #mu_est = max(min(100.0 * norm(get_primal_res(iter),Inf), iter.point.mu, mean(iter.point.y .* iter.point.s) ), feas_obj) + norm(get_primal_res(iter),Inf) * 0.01
               #mu_est = mean(iter.point.y .* iter.point.s)
               #mu_est = (mean(iter.point.y .* iter.point.s) + mean(-iter.point.y .* get_primal_res(iter))) / 4.0
               mu_ub = (mean(iter.point.y) + 1.0) * mean(-get_primal_res(iter))
               #dual_avg = norm(eval_grad_lag(iter),1) / length(iter.point.mu)
               dual_avg = scaled_dual_feas(iter, pars)
               dual_progress = dual_avg < norm(get_primal_res(iter), Inf)
               # * 10.0
               #dual_progress = dual_avg < mu_ub * 10.0
               delta_small = get_delta(iter) < get_mu(iter) * (1.0 + norm(get_y(iter),Inf))
               lag_grad = norm(eval_grad_lag(iter),Inf) < norm(get_grad(iter),Inf) + (norm(get_primal_res(iter), Inf) + 1.0) #+ sqrt(norm(get_y(iter),Inf))

               be_aggressive = is_feas && (lag_grad || !pars.lag_grad_test) && dual_progress && (delta_small || !pars.inertia_test)


               if be_aggressive
                     #=prm_tol_ratio = pars.tol / norm(get_primal_res(iter),Inf)
                     if prm_tol_ratio >= 1.0 && false
                        step_type = "mu"
                        reduct_factors = Class_reduction_factors(0.8, 0.0, 1e-1)
                     elseif false  #10.0 * mu_P_ratio^(7/4)
                       step_type = "prm"
                       reduct_factors = Class_reduction_factors(1e-1, 0.0, 0.8)
                     else
                     end=#
                     if norm(get_primal_res(iter),Inf) > pars.tol || !pars.pause_primal
                       step_type = "agg"
                       reduct_factors = Reduct_affine()
                     else
                       reduct_factors = Class_reduction_factors(1.0, 0.0, 0.0)
                       step_type = "mu"
                     end
                     #reduct_factors = Class_reduction_factors(0.5, 0.5, 0.5)
                     ls_mode = pars.ls_mode_agg;

                     q = pars.min_step_size_agg_ratio * min(1.0, 1.0 / maximum(- get_primal_res(iter) ./ iter.point.s))
                     actual_min_step_size = q
                     #@show q
                     iter.point.mu = mu_est
                     centre_dual!(iter.point, pars.comp_feas_agg)
                 else
                      reduct_factors = Reduct_stable()
                      step_type = "stb"
                      ls_mode = pars.ls_mode_stable_correction;
                      actual_min_step_size = pars.min_step_size_stable
                 end

               if i == 1
                   update_H!(iter, timer)

                   ipopt_strategy!(iter, kkt_solver, pars, timer)

                   start_advanced_timer(timer, "STEP/first")

                   for i = 1:100
                     success, new_iter, ls_info = take_step!(iter, reduct_factors, kkt_solver, ls_mode, filter, pars, actual_min_step_size, timer)

                     if success == :success
                       iter = new_iter
                       break
                     elseif i < 100
                       set_delta(iter, max(get_delta(iter) * 8.0, 1e-8))
                       factor!(kkt_solver, get_delta(iter), timer)
                     else
                       pause_advanced_timer(timer, "STEP/first")
                       println("Terminated due to max delta")
                       println("delta=$(get_delta(iter)), step_type=$step_type, min_step_size=$actual_min_step_size")
                       @show ls_info
                       return iter, :MAX_DELTA, progress, t, false
                     end
                   end

                   pause_advanced_timer(timer, "STEP/first")
               else
                     if be_aggressive
                       iter.point.mu = mu_est
                       centre_dual!(iter.point, pars.comp_feas_agg)
                     end

                     start_advanced_timer(timer, "STEP/correction")
                     success, new_iter, ls_info = take_step!(iter, reduct_factors, kkt_solver, ls_mode, filter, pars, actual_min_step_size, timer)

                     if success == :success
                       iter = new_iter
                       if step_type == "agg"
                         dir_size_agg = norm(kkt_solver.dir.x, 2)
                         update_prox!(iter, pars)
                       end
                     end

                     pause_advanced_timer(timer, "STEP/correction")
               end

               add!(filter, iter, pars)
               record_progress!(t,  step_type, iter, kkt_solver, ls_info, reduct_factors, progress, pars)

               if success != :success
                 break
               end
             end
      end

    #=catch(e)
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
