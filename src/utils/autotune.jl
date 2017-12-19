function autotune(nlp_raw, pars::Class_parameters)
    nlp = Class_CUTEst(nlp_raw);

    best_it = Inf
    best_mu_ratio = 0.0
    success = true

    for i = -2:5
      local_timer = class_advanced_timer()
      start_advanced_timer(local_timer)
      mu_ratio = 10.0^(i)
      pars.init.mu_scale = mu_ratio
      println("--------  MU_SCALE = $(mu_ratio) ----------")
      start_advanced_timer(local_timer, "INIT");
      intial_it = init(nlp, pars, local_timer);
      pause_advanced_timer(local_timer, "INIT");

      @assert(is_feasible(intial_it, pars.ls.comp_feas))
      iter, status, hist, t, err = one_phase_IPM(intial_it, pars, local_timer);
      if t < best_it - 1 && status == :optimal
        best_it = t
        best_mu_ratio = mu_ratio
      end

      pause_advanced_timer(local_timer)
      print_timer_stats(local_timer)
    end

    return success, best_it, best_mu_ratio
end
