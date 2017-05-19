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
