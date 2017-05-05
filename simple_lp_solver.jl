function simple_LP_solver(intial_it::Class_iterate, pars::Class_parameters)
  start_advanced_timer()

  iter = intial_it;

  kkt_solver = pick_KKT_solver(pars);
  initialize!(kkt_solver, intial_it)

  progress = Array{alg_history,1}()

  head_progress()
  record_progress_first_it!(iter, kkt_solver, progress, pars)


  for t = 1:300
     factor!(kkt_solver, iter)

     if scaled_dual_feas(iter) < get_mu(iter)
       step_type = "agg"
       reduct_factors = Reduct_affine()
       # get_mu(it) * (1.0 - eta)

       tmp_iter, ls_info = take_step!(iter, reduct_factors, kkt_solver, pars)
       step_size = ls_info.step_size_P
       eta = min(0.99, max((1.0 - step_size)^2, 0.001))

       reduct_factors = Eta_reduct(eta, :symmetric)
       iter, ls_info = take_step!(iter, reduct_factors, kkt_solver, pars)
     else
        step_type = "sta"
        reduct_factors = Reduct_stable()
        iter, ls_info = take_step!(iter, reduct_factors, kkt_solver, pars)
     end

     record_progress!(t, step_type, iter, kkt_solver, ls_info, reduct_factors, progress, pars)


     status = terminate(iter, pars)
     if status != false
       println("Terminated with ", status)
       return iter, status, progress
     end
  end
end

function better_LP_solver(intial_it::Class_iterate, pars::Class_parameters)
  start_advanced_timer()

  iter = intial_it;

  kkt_solver = pick_KKT_solver(pars);
  initialize!(kkt_solver, intial_it)

  progress = Array{alg_history,1}()

  head_progress()
  record_progress_first_it!(iter, kkt_solver, progress, pars)
  init_step_size = 1.0

  success = :success

  for t = 1:pars.max_it
     if init_step_size < 1e-2
       set_delta(iter, get_delta(iter) * 10.0)
     end
     factor!(kkt_solver, iter)

     for i = 1:3
       if i > 1 && scaled_dual_feas(iter) < get_mu(iter) && norm(comp(iter),Inf) < 0.99 * get_mu(iter)
         step_type = "agg"
         reduct_factors = Reduct_affine()
         accept = accept_aggressive;
       elseif true
         step_type = "stb"
         reduct_factors = Reduct_stable()
         accept = accept_aggressive;
         #accept = accept_stable;
       else
         break
       end

       success, iter, ls_info = take_step!(iter, reduct_factors, kkt_solver, accept, pars)

       record_progress!(t, step_type, iter, kkt_solver, ls_info, reduct_factors, progress, pars)

       if i == 1
          init_step_size = ls_info.step_size_P
          if success != :success
            set_delta(iter, max(1e-8,get_delta(iter) * 10.0))
          end
       end

       if ls_info.step_size_P < 1e-1 * init_step_size || success != :success
         break
       end

       status = terminate(iter, pars)
       if status != false
         println("Terminated with ", status)
         return iter, status, progress
       end
     end
  end

  println("Max iterations reached")
  return iter, :max_it, progress
end
