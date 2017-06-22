

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
      new_iter = nothing;
      agg_next_step = false;
      dir_size_agg = Inf
      dir_size_stable = 0.0

      start_time = time()

      for t = 1:pars.max_it
            #@show norm(get_s(iter),2)
            #@show norm(get_y(iter),2)
            #@show norm(get_y(iter),Inf)
            #@show maximum(a_neg_part .* y)
             #g = eval_grad_lag(nlp, iter.point.x, zeros(ncon(iter)))
             #@show mean(abs(g))
             @assert(is_feasible(iter, pars.comp_feas))

             #@show

             for i = 1:(pars.max_it_corrections+1)
               #iter.cache.J = eval_jac(iter.nlp, iter.point.x)
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

               if pars.threshold_type == :mu
                 threshold = pars.aggressive_dual_threshold * get_mu(iter)
               elseif pars.threshold_type == :mu_primal
                  threshold = pars.aggressive_dual_threshold * min(get_mu(iter),norm(get_primal_res(iter),Inf))
               elseif pars.threshold_type == :primal
                   threshold = pars.aggressive_dual_threshold * norm(get_primal_res(iter),Inf)
               else
                 error("threshold_type not known")
               end

               is_feas = is_feasible(iter, pars.comp_feas_agg)
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

               #* (1.0 + norm(iter.point.x,2))
               #dual_progress = norm(eval_grad_lag(iter),2) < dot(iter.point.s - get_primal_res(iter), iter.point.y)
               #dual_progress = dual_avg < mu_est #+ sqrt(norm(get_primal_res(iter),Inf))
               dual_progress = dual_avg < norm(get_primal_res(iter), Inf) * 10.0
               #dual_progress = dual_avg < mu_ub * 10.0
               delta_small = get_delta(iter) < sqrt(get_mu(iter)) * (1.0 + norm(get_y(iter),Inf))
               lag_grad = norm(eval_grad_lag(iter),Inf) < norm(get_grad(iter),Inf) + (norm(get_primal_res(iter), Inf) + 1.0) #+ sqrt(norm(get_y(iter),Inf))

               be_aggressive = is_feas && (lag_grad || !pars.lag_grad_test) && dual_progress && (delta_small || !pars.inertia_test)
               #@show norm(get_grad(iter),Inf)
               #=if false
                   else
                  mu = get_mu(iter)
                  be_aggressive = is_feas && scaled_dual_feas(iter, pars) < mu && norm(eval_grad_lag(iter),Inf) < norm(get_grad(iter),Inf) / 2 && (delta_small || !pars.inertia_test)
               end=#

               if i == 1
                   step_type = "stb"
                   reduct_factors = Reduct_stable()

                   update_H!(iter, timer)

                   start_advanced_timer(timer, "STEP/first-stable")
                   #success, new_iter, ls_info = optimal_stable_step!(iter, reduct_factors, kkt_solver, pars.ls_mode_stable, filter, pars)
                   if pars.trust_region
                     success, new_iter, ls_info = trust_region_strategy!(iter,kkt_solver, filter, pars, timer)
                   else
                     success, new_iter, ls_info = ipopt_strategy!(iter, reduct_factors, kkt_solver, filter, pars, timer)
                   end
                   pause_advanced_timer(timer, "STEP/first-stable")


                   init_step_size = ls_info.step_size_P
                   dir_size_stable = norm(kkt_solver.dir.x, 2)

                   if pars.eigen_search && ls_info.step_size_P == 1.0 && get_delta(iter) > norm(comp(new_iter),Inf) + scaled_dual_feas(new_iter, pars)
                      min_eig = compute_eigenvector!(kkt_solver, new_iter, timer)
                      @show min_eig
                      dir = scale_direction(kkt_solver.dir, dir_size_stable / 10.0)
                      status, new_iter = eigenvector_ls(new_iter, dir, pars)
                   end

                   if pars.stb_before_agg || !be_aggressive
                     iter = new_iter
                   end
               else

                   #ratio = minimum(iter.point.s ./ (iter.point.primal_scale * iter.primal_residual_intial))
                   #@show ratio


                   if be_aggressive #&& i == 2
                     #=g = eval_grad_lag(iter.nlp, iter.point.x, zeros(length(iter.point.y)))

                     target_mu_reduction = (norm(g, Inf) + 1e-2) * norm(eval_primal_residual(iter),Inf) / get_mu(iter)
                     if target_mu_reduction > 1e1
                       step_type = "prm"
                       reduct_factors = Class_reduction_factors(1.0 / target_mu_reduction, 0.0, 1.0)
                     elseif target_mu_reduction < 1e-1
                       step_type = "mu"
                       reduct_factors = Class_reduction_factors(1.0, 0.0, target_mu_reduction)
                     else
                         step_type = "agg"
                         reduct_factors = Reduct_affine()
                      end=#

                       #reduct_factors = Reduct_affine()
                    #mode = :no_change

                    #norm(get_grad(iter),Inf) /

                    if true
                     prm_tol_ratio = pars.tol / norm(get_primal_res(iter),Inf)
                     if prm_tol_ratio >= 1.0 && false
                        step_type = "mu"
                        reduct_factors = Class_reduction_factors(0.8, 0.0, 1e-1)
                     elseif false  #10.0 * mu_P_ratio^(7/4)
                       step_type = "prm"
                       reduct_factors = Class_reduction_factors(1e-1, 0.0, 0.8)
                     else
                        step_type = "agg"
                        reduct_factors = Reduct_affine() #
                        #reduct_factors = Class_reduction_factors(prm_tol_ratio, 0.0, 0.0 )
                     end
                   end

                    mu_P_ratio = iter.point.mu / primal_inf
                    seems_to_be_infeasible = false
                    feas_with_big_mu = false

                    if pars.adaptive_mu == :test1
                       #
                       prm_to_dl = (1.0 + mean(iter.point.y)) / (mean(iter.point.s) + 1.0) #(norm(iter.point.y, Inf) + 1.0) / (norm(iter.point.x, Inf) + 1.0)
                       seems_to_be_infeasible = prm_to_dl > 20.0 * mu_P_ratio^(3/2) # 10.0 * ((5.0 + norm(g,Inf)) *
                       feas_with_big_mu = prm_to_dl < 0.1 * mu_P_ratio^(1/2)
                     elseif pars.adaptive_mu == :test2
                       feas = norm(iter.point.primal_scale * iter.primal_residual_intial, Inf)
                       #g = eval_grad_lag(iter.nlp, iter.point.x, zeros(length(iter.point.y)))
                       mu = iter.point.mu
                       p_scale = iter.point.primal_scale
                       seems_to_be_infeasible = norm(eval_farkas(iter),Inf) < 10.0 * p_scale^(3/2)
                       feas_with_big_mu = false  #norm(eval_farkas(iter),Inf) > 0.1 * p_scale^(3/2) && norm(iter.point.y, Inf) > 1.0 / min(0.01, sqrt(mu))
                     elseif pars.adaptive_mu == :test3
                       #prm_to_dl = mean(iter.point.y) + mean(iter.point.s) #norm(iter.point.y, Inf) + norm(iter.point.s, Inf) #+ 1.0) / (norm(iter.point.s, Inf) + 1.0)
                       prm_to_dl = mean(iter.point.y) # mean([iter.point.y; iter.point.s])
                       #prm_to_dl = median(iter.point.y) / median(iter.point.s)
                       #prm_to_dl = median(iter.point.y) / 2 + median(iter.point.s) / 2
                       @show norm(get_grad(iter),Inf)
                       seems_to_be_infeasible = prm_to_dl > 100.0 * mu_P_ratio^(5/4) # 10.0 * ((5.0 + norm(g,Inf)) *
                       feas_with_big_mu = prm_to_dl < 0.1 * (mu_P_ratio)^(3/4)
                     elseif pars.adaptive_mu == :test4
                       seems_to_be_infeasible = iter.point.mu < 0.1 * mu_est # 10.0 * ((5.0 + norm(g,Inf)) *
                       feas_with_big_mu = iter.point.mu > 10.0 * mu_est
                     end


                       if be_aggressive
                         #iter.point.mu = dot(iter.point.y, iter.point.s - get_primal_res(iter)) / (2.0 * length(iter.point.y))
                         #@show mean(iter.point.y .* -get_primal_res(iter)), mean(iter.point.y) * mean(abs(get_primal_res(iter)))

                         iter.point.mu = mu_est
                         #mean(iter.point.y) * mean(abs(get_primal_res(iter))) / 10.0
                         centre_dual!(iter.point, pars.comp_feas_agg)

                         #record_progress!(t,  "Δmu", iter, kkt_solver, ls_info, reduct_factors, progress, pars)
                       end

                       #step_type = "agg"
                       #reduct_factors = Reduct_affine()

                       #reduct_factors = Class_reduction_factors(0.2, 0.0, 0.0)
                       ls_mode = pars.ls_mode_agg;
                       actual_min_step_size = 0.0
                   elseif true
                     #if i == 2
                     #  iter = new_iter
                     #end
                     step_type = "stb"
                     reduct_factors = Reduct_stable()
                     #accept = accept_aggressive;
                     ls_mode = pars.ls_mode_stable_correction;
                     actual_min_step_size = pars.min_step_size_stable
                   else
                     break
                   end

                   start_advanced_timer(timer, "STEP/correction")
                   success, iter, ls_info = take_step!(iter, reduct_factors, kkt_solver, ls_mode, filter, pars, actual_min_step_size, timer)
                   if step_type == "agg"
                     dir_size_agg = norm(kkt_solver.dir.x, 2)
                     update_prox!(iter, pars)


                    #println("mu_est = ",)
                    #mu_est =  mean(iter.point.s .* iter.point.y)
                    #iter.point.y = iter.point.y * (mu_est / iter.point.mu)
                    #iter.point.mu = mu_est
                    if false
                      if feas_with_big_mu #|| norm(get_primal_res(iter),Inf) < pars.tol
                          record_progress!(t,  "dmu", iter, kkt_solver, ls_info, reduct_factors, progress, pars)
                          iter.point.y *= 1/4.0
                          iter.point.mu *= 1/4.0
                       elseif seems_to_be_infeasible  #10.0 * mu_P_ratio^(7/4)

                         #println("increase mu")
                         iter.point.y *= 4.0
                         iter.point.mu *= 4.0

                         record_progress!(t,  "imu", iter, kkt_solver, ls_info, reduct_factors, progress, pars)
                       end
                     end
                   end

                   pause_advanced_timer(timer, "STEP/correction")
               end

               add!(filter, iter, pars)
               record_progress!(t,  step_type, iter, kkt_solver, ls_info, reduct_factors, progress, pars)
               #g = eval_grad_lag(iter.nlp, iter.point.x, zeros(length(iter.point.y)))
               #J = eval_jac(iter)
               #@show norm(J, Inf)
               #H = eval_lag_hess(iter)
               #@show norm(H, Inf)



               #=y_av = mean(iter.point.y)
               s_av = mean(iter.point.s)
               ratios = min(iter.point.y / y_av, iter.point.s / s_av) / iter.point.mu
               @show mean(ratios)
               @show minimum(ratios)
               @show maximum(ratios)=#
               #@show norm(g, Inf)
               #ratio = mean(iter.point.y + iter.point.s) * (iter.point.mu / min(minimum(iter.point.y),minimum(iter.point.s)))
               #@show ratio

               #@show maximum(1.0 ./ (iter.point.y + iter.point.s))
               #=
               Dμ = abs((kkt_solver.dir.s ./ iter.point.s) .* (kkt_solver.dir.y ./ iter.point.y))
               ratio_thing = maximum(Dμ / get_mu(iter))

               @show ratio_thing
               @show norm(eval_grad_lag(iter.nlp,iter.point.x, zeros(length(iter.point.y))),Inf)
               if ratio_thing < 10.0
                 agg_next_step = true
               else
                 agg_next_step = false
               end=#

               #@show strict_comp = maximum( (iter.point.s .* iter.point.y) ./ min(iter.point.s, iter.point.y) )
               #if iter.point.mu < 1.0 && strict_comp > 100.0
               #  my_warn("strict complemenarity does not seem to hold")
               #end

               if i > 1 && ls_info.step_size_P < pars.min_step_size_correction || success != :success
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
