

function head_progress()
  println( pd("it",5), pd("step", 5), pd("eta"), pd("α_P"), pd("α_D"), pd("ls",3), pd("|dx|"), pd("|dy|"), pd("kkt"), "|", pd("mu"), pd("dual"), pd("primal"),pd("comp"), pd("inf"), "|", pd("delta"), pd("|y|"), pd("|x|"), pd("∇phi"), pd("phi") )
end

type alg_history
    t::Int64
    step_type::String
    eta::Class_reduction_factors
    ls_info::abstract_ls_info
    dir_x_norm::Float64
    dir_y_norm::Float64
    dir_s_norm::Float64
    kkt_err::Class_kkt_error
    mu::Float64
    norm_grad_lag::Float64
    primal_residual::Float64
    comp::Float64
    dot_sy::Float64
    farkas::Float64
    delta::Float64
    eval_merit_function::Float64
    eval_phi::Float64
    eval_grad_phi::Float64
    y_norm::Float64
    x_norm::Float64

    function alg_history()
      return new()
    end
end

function record_progress_first_it!(iter::Class_iterate, kss::abstract_KKT_system_solver, progress::Array{alg_history,1}, pars::Class_parameters)
    record_progress!(0, "it0", iter, kss, Blank_ls_info(), Class_reduction_factors(), progress, pars)
end


function record_progress!(t::Int64, step_type::String, iter::Class_iterate, kss::abstract_KKT_system_solver, ls_info::abstract_ls_info, eta::Class_reduction_factors, progress::Array{alg_history,1}, par::Class_parameters)
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
    hist.norm_grad_lag = norm(eval_grad_lag(iter),Inf)
    hist.primal_residual = norm(eval_primal_residual(iter),Inf)
    hist.comp = norm(comp(iter),Inf)
    hist.farkas = norm(eval_farkas(iter),Inf)
    hist.delta = get_delta(iter)
    hist.eval_merit_function = eval_merit_function(iter, par)
    hist.eval_phi = eval_phi(iter)
    hist.eval_grad_phi = norm(eval_grad_phi(iter),Inf)
    hist.y_norm = norm(get_y(iter),Inf)
    hist.x_norm = norm(get_x(iter),Inf)

    #@show hist.eta.P, hist.eta.mu, hist.eta.P

    kkt_ratio = hist.kkt_err.ratio

    push!(progress, hist)

    if par.output_level >= 2
      println(pd(t,5), pd(step_type, 5), rd(eta.mu),  rd(ls_info.step_size_P),  rd(ls_info.step_size_D), pd(ls_info.num_steps,3), rd(hist.dir_x_norm), rd(hist.dir_y_norm), rd(kkt_ratio), "|", rd(hist.mu), rd(hist.norm_grad_lag), rd(hist.primal_residual), rd(hist.comp), rd(hist.farkas), "|", rd(hist.delta),  rd(hist.y_norm), rd(hist.x_norm), rd(hist.eval_grad_phi), rd(hist.eval_merit_function, 5))
    end
end


function Reduct_affine()
  return Class_reduction_factors(0.0, 0.0, 0.0);
end

function Reduct_stable()
  return Class_reduction_factors(1.0, 0.0, 1.0);
end

function Eta_reduct(eta::Float64, strategy::Symbol)
    if strategy == :dual_agg
      return Class_reduction_factors(eta, 0.0, eta);
    elseif strategy == :symmetric
      return Class_reduction_factors(eta, eta, eta);
    else
      error("This eta reduction strategy does not exist")
    end
end

function dual_scale(iter::Class_iterate)
    return 100.0 / max(norm(get_y(iter), Inf),100.0)
end

function scaled_dual_feas(iter::Class_iterate)
    return norm(eval_grad_lag(iter),Inf) * dual_scale(iter)
end

function is_feasible(it::Class_iterate, comp_feas::Float64)
    mu = get_mu(it);
    Sy = get_s(it) .* get_y(it)
    if length(Sy) > 0
      return all(it.point.s .> 0.0) && all(it.point.y .> 0.0) && maximum(Sy) / mu  <= 1.0 / comp_feas && minimum(Sy) / mu >= comp_feas
    else
      return true
    end
end

function terminate(iter::Class_iterate, par::Class_parameters)
    tol = par.tol

    if scaled_dual_feas(iter) < tol && get_mu(iter) < tol && norm(eval_primal_residual(iter),Inf) < tol
        return :optimal
    elseif norm(eval_farkas(iter),Inf) < tol && norm(eval_primal_residual(iter),Inf) > tol && norm(get_y(iter),Inf) > 1/tol
        return :primal_infeasible
    elseif eval_f(iter) < -1.0/tol^2 # TMP!!!
        return :dual_infeasible
    else
        return false
    end
end

function take_step!(iter::Class_iterate, eta::Class_reduction_factors, kkt_solver::abstract_KKT_system_solver, ls_mode::Symbol, filter::Array{Class_filter,1}, pars::Class_parameters, min_step_size::Float64)
    kkt_associate_rhs!(kkt_solver, iter, eta)
    compute_direction!(kkt_solver)

    if false #kkt_solver.kkt_err_norm.ratio > 1e-1
        @show kkt_solver.kkt_err_norm
        return :failure, iter, Class_ls_info()
    end

    success, iter, step_info = simple_ls(iter, kkt_solver.dir, ls_mode, filter, pars, min_step_size)

    return success, iter, step_info
end
