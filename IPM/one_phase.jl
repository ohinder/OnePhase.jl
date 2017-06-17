

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

               #inertia_good = get_delta(iter) < max(get_mu(iter), 1.0)
               primal_inf = norm(iter.primal_residual_intial, Inf) * iter.point.primal_scale
               #min(primal_inf, get_mu(iter))
               threshold = pars.aggressive_dual_threshold * min(get_mu(iter),norm(get_primal_res(iter),Inf))
               dual_progress = scaled_dual_feas(iter, pars) < threshold && is_feasible(iter, pars.comp_feas_agg)
               delta_small = get_delta(iter) < sqrt(get_mu(iter)) * (1.0 + norm(get_y(iter),Inf))
               be_aggressive = dual_progress && (delta_small || !pars.inertia_test)

               if i == 1
                   step_type = "stb"
                   reduct_factors = Reduct_stable()

                   update_H!(iter, timer)

                   start_advanced_timer(timer, "STEP/first-stable")
                   #success, new_iter, ls_info = optimal_stable_step!(iter, reduct_factors, kkt_solver, pars.ls_mode_stable, filter, pars)
                   success, new_iter, ls_info = ipopt_strategy!(iter, reduct_factors, kkt_solver, filter, pars, timer)
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
                    if pars.adaptive_mu == :test1
                       mu_P_ratio = iter.point.mu / primal_inf
                       #
                       seems_to_be_infeasible = mu_P_ratio < 1e5 && norm(iter.point.y, Inf) > 10.0 * mu_P_ratio^(3/2) # 10.0 * ((5.0 + norm(g,Inf)) *
                       feas_with_big_mu = norm(iter.point.y, Inf) < 0.1 * mu_P_ratio^(1/2)
                     elseif pars.adaptive_mu == :test2
                       feas = norm(iter.point.primal_scale * iter.primal_residual_intial, Inf)
                       #g = eval_grad_lag(iter.nlp, iter.point.x, zeros(length(iter.point.y)))
                       mu = iter.point.mu
                       p_scale = iter.point.primal_scale
                       seems_to_be_infeasible = norm(eval_farkas(iter),Inf) < 10.0 * p_scale^(3/2)
                       feas_with_big_mu = false  #norm(eval_farkas(iter),Inf) > 0.1 * p_scale^(3/2) && norm(iter.point.y, Inf) > 1.0 / min(0.01, sqrt(mu))
                     elseif pars.adaptive_mu == :none
                       seems_to_be_infeasible = false
                       feas_with_big_mu = false
                     end


                      if true
                       if norm(get_primal_res(iter),Inf) < pars.tol
                          step_type = "mu"
                          reduct_factors = Class_reduction_factors(0.8, 0.0, 1e-1)
                       elseif false  #10.0 * mu_P_ratio^(7/4)
                         step_type = "prm"
                         reduct_factors = Class_reduction_factors(1e-1, 0.0, 0.8)
                       else
                          step_type = "agg"
                          #reduct_factors = Reduct_affine() #
                          reduct_factors = Class_reduction_factors(0.5, 0.5, 0.5) #Reduct_affine()
                       end
                     end

                     if feas_with_big_mu #|| norm(get_primal_res(iter),Inf) < pars.tol
                        iter.point.y *= 0.5
                        iter.point.mu *= 0.5
                     elseif seems_to_be_infeasible  #10.0 * mu_P_ratio^(7/4)
                       iter.point.y *= 2.0
                       iter.point.mu *= 2.0
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
