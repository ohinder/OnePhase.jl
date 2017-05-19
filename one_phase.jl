function one_phase_IPM(intial_it::Class_iterate, pars::Class_parameters)
  start_advanced_timer()

  iter = intial_it;

  filter = Array{Class_filter,1}();
  kkt_solver = pick_KKT_solver(pars);
  initialize!(kkt_solver, intial_it)

  progress = Array{alg_history,1}()

  head_progress()
  record_progress_first_it!(iter, kkt_solver, progress, pars)
  init_step_size = 1.0

  success = :success
  new_iter = nothing;

  for t = 1:pars.max_it
     for i = 1:(pars.max_it_corrections+1)
       if i == 1
           step_type = "stb"
           reduct_factors = Reduct_stable()
           #success, iter, ls_info = optimal_stable_step!(iter, reduct_factors, kkt_solver, pars)
           success, iter, ls_info = ipopt_strategy!(iter, reduct_factors, kkt_solver, pars.ls_mode_stable, filter, pars)
           init_step_size = ls_info.step_size_P
       else
           #inertia_good = get_delta(iter) < max(get_mu(iter), 1.0)
           dual_progress = scaled_dual_feas(iter) < get_mu(iter) && is_feasible(iter,pars.comp_feas_agg)
           be_aggressive = dual_progress

           if be_aggressive #&& i == 2
             step_type = "agg"
             reduct_factors = Reduct_affine()
             ls_mode = pars.ls_mode_agg;
             actual_min_step_size = 0.0
           elseif true
             #if i == 2
             #  iter = new_iter
             #end
             step_type = "stb"
             reduct_factors = Reduct_stable()
             #accept = accept_aggressive;
             ls_mode = pars.ls_mode_stable;
             actual_min_step_size = pars.min_step_size_stable
           else
             break
           end

           success, iter, ls_info = take_step!(iter, reduct_factors, kkt_solver, ls_mode, filter, pars, actual_min_step_size)
       end

       add!(filter, iter)
       record_progress!(t,  step_type, iter, kkt_solver, ls_info, reduct_factors, progress, pars)
       #@show maximum(1.0 ./ (iter.point.y + iter.point.s))

       if ls_info.step_size_P < pars.min_step_size_correction || success != :success
         break
       end

       status = terminate(iter, pars)
       if status != false
         println("Terminated with ", status)
         return iter, status, progress, t
       end
     end
  end

  println("Max iterations reached")
  return iter, :max_it, progress, t
end
