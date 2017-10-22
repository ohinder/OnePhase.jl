type alg_history <: abstract_alg_history
    t::Int64
    step_type::String
    eta::Class_reduction_factors
    ls_info::abstract_ls_info
    dir_x_norm::Float64
    dir_y_norm::Float64
    dir_s_norm::Float64
    kkt_err::Class_kkt_error
    mu::Float64
    fval::Float64
    dual_scaled::Float64
    norm_grad_lag::Float64
    primal_residual::Float64
    con_vio::Float64
    comp::Float64
    comp_ratio::Float64
    dot_sy::Float64
    farkas::Float64
    delta::Float64
    eval_merit_function::Float64
    eval_phi::Float64
    eval_grad_phi::Float64
    y_norm::Float64
    x_norm::Float64
    tot_num_fac::Int64
    num_fac_inertia::Int64

    function alg_history()
      return new()
    end
end

type ipopt_alg_history <: abstract_alg_history
  t::Int64
  #mu::Float64
  fval::Float64
  #dual_scaled::Float64
  norm_grad_lag::Float64
  comp::Float64
  #comp_ratio::Float64
  #dot_sy::Float64
  #primal_residual::Float64
  #farkas::Float64
  #delta::Float64
  #eval_merit_function::Float64
  #eval_phi::Float64
  #eval_grad_phi::Float64
  con_vio::Float64
  y_norm::Float64
  x_norm::Float64

  function ipopt_alg_history()
    return new()
  end
end

function head_progress()
  println( pd("it",3), pd("step", 6), pd("eta",8), pd("α_P"), pd("α_D"), pd("ls",3), pd("|dx|"), pd("|dy|"), pd("N err"), "|", pd("mu"), pd("dual"), pd("primal"), pd("cmp scaled"), pd("inf"), "|", pd("delta"), pd("#ch",3), pd("|x|"), pd("|y|"), pd("∇phi"), pd("phi") )
end

function record_progress_first_it!(iter::Class_iterate, kss::abstract_KKT_system_solver, progress::Array{alg_history,1}, pars::Class_parameters)
    record_progress!(0, "it0", iter, kss, Blank_ls_info(), Class_reduction_factors(), progress, 0, 0, pars)
end


function record_progress!(t::Int64, step_type::String, iter::Class_iterate, kss::abstract_KKT_system_solver, ls_info::abstract_ls_info, eta::Class_reduction_factors, progress::Array{alg_history,1}, num_fac_inertia::Int64, tot_num_fac::Int64, par::Class_parameters)
    hist = alg_history()
    hist.t = t;
    hist.step_type = step_type;
    hist.eta = eta;
    hist.ls_info = ls_info;
    hist.dir_x_norm = norm(kss.dir.x,Inf);
    hist.dir_y_norm = norm(kss.dir.y,Inf);
    hist.dir_s_norm = norm(kss.dir.s,Inf);
    hist.kkt_err = kss.kkt_err_norm
    hist.mu = get_mu(iter)
    hist.dot_sy = dot(iter.point.s,iter.point.y)/length(iter.point.s)
    hist.norm_grad_lag = norm(eval_grad_lag(iter, iter.point.mu),Inf)
    hist.dual_scaled = scaled_dual_feas(iter, par)
    hist.primal_residual = norm(get_primal_res(iter),Inf)
    hist.con_vio = get_max_vio(iter)
    hist.comp = norm(comp(iter),Inf)
    hist.comp_ratio = comp_ratio_max(iter) * par.comp_feas_agg #norm(comp(iter),Inf)
    hist.farkas = eval_farkas(iter)
    hist.delta = get_delta(iter)
    hist.fval = get_fval(iter)
    hist.eval_merit_function = eval_merit_function(iter, par)
    hist.eval_phi = eval_phi(iter,iter.point.mu)
    hist.eval_grad_phi = norm(eval_grad_phi(iter, iter.point.mu),Inf)
    hist.y_norm = norm(get_y(iter),Inf)
    hist.x_norm = norm(get_x(iter),Inf)
    hist.tot_num_fac = tot_num_fac
    hist.num_fac_inertia = num_fac_inertia

    #@show hist.eta.P, hist.eta.mu, hist.eta.P

    kkt_ratio = hist.kkt_err.ratio

    push!(progress, hist)

    if par.output_level >= 2
      println(pd(t,5), pd(step_type[1], 3), rd(eta.mu),  rd(ls_info.step_size_P),  rd(ls_info.step_size_D), pd(ls_info.num_steps,2), rd(hist.dir_x_norm), rd(hist.dir_y_norm), rd(kkt_ratio), "|", rd(hist.mu), rd(hist.dual_scaled), rd(hist.primal_residual), rd(hist.comp_ratio), rd(hist.farkas), "|", rd(hist.delta), pd(hist.tot_num_fac,3), rd(hist.x_norm), rd(hist.y_norm), rd(hist.eval_grad_phi), rd(hist.eval_merit_function, 5))
    end
end
