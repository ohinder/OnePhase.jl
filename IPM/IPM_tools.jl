function centre_dual!(point::Class_point, comp_feas::Float64)
  y_c = point.mu ./ point.s

  point.y = min( y_c / comp_feas, max(point.y, comp_feas * y_c))
end


function dual_scale(iter::Class_iterate, pars::Class_parameters)
    if pars.dual_scale_mode == :scaled
      return pars.dual_scale_threshold / max(norm(get_y(iter), Inf), pars.dual_scale_threshold)
    elseif pars.dual_scale_mode == :sqrt
      return pars.dual_scale_threshold / max(sqrt(norm(get_y(iter), Inf)), pars.dual_scale_threshold)
    elseif pars.dual_scale_mode == :exact
      return 1.0
    elseif pars.dual_scale_mode == :primal_dual
      return pars.dual_scale_threshold / max(sqrt(norm(get_y(iter), Inf) * norm(get_s(iter), Inf)), + pars.dual_scale_threshold)
    else
      throw("dual scale type does not exist")
    end
end

function scaled_dual_feas(iter::Class_iterate, pars::Class_parameters)
    return norm(eval_grad_lag(iter),Inf) * dual_scale(iter, pars)
end

function check_for_nan(point::Class_point)
  if isbad(point.s)
    error("NaN in s")
  end

  if isbad(point.y)
    #@show point.y, point.s
    error("NaN in y")
  end

  if isbad(point.x)
    error("NaN in x")
  end

  if isbad(point.mu)
    error("NaN in mu")
  end
end

function is_feasible(it::Class_iterate, comp_feas::Float64)
    check_for_nan(it.point)

    mu = get_mu(it);
    s =  get_s(it)
    y =  get_y(it)

    Sy = s .* y
    if length(Sy) > 0
      return all(s .> 0.0) && all(y .> 0.0) && maximum(Sy) / mu  <= 1.0 / comp_feas && minimum(Sy) / mu >= comp_feas
    else
      return true
    end
end

function terminate(iter::Class_iterate, par::Class_parameters)
    tol = par.tol
    #a_neg_part = max(-eval_a(iter),0.0)
    #y = get_y(iter)
    #J = eval_jac(iter)
    fark_feas1 = eval_farkas(iter)
    #fark_feas1 = norm(y' * J, Inf) / min(norm(y, Inf), maximum(a_neg_part .* y))
    #fark_feas2 = norm(eval_grad_lag(iter, 0.0), Inf) /  min(norm(y, Inf), maximum(a_neg_part .* y))
    #fark_feas3 = norm(eval_grad_lag(iter, 1.0), Inf) /  min(norm(y, Inf), maximum(a_neg_part .* y))

    #@show fark_feas1, fark_feas2, fark_feas3
    max_vio = get_max_vio(iter)

    comp_scaled = maximum(iter.point.s .* iter.point.y) * dual_scale(iter, par)

    if scaled_dual_feas(iter, par) < tol && comp_scaled < tol && max_vio < tol
        return :optimal
    elseif fark_feas1 < tol && max_vio > tol && norm(get_y(iter), Inf) > 1/tol #&& norm(get_y(iter),Inf) > 1/tol
        return :primal_infeasible
    elseif get_fval(iter) < -1.0/tol^2 # TMP!!!
        return :dual_infeasible
    else
        return false
    end
end

function take_step!(iter::Class_iterate, eta::Class_reduction_factors, kkt_solver::abstract_KKT_system_solver, ls_mode::Symbol, filter::Array{Class_filter,1}, pars::Class_parameters, min_step_size::Float64, timer::class_advanced_timer)
    kkt_associate_rhs!(kkt_solver, iter, eta, timer)
    compute_direction!(kkt_solver, timer)

    if false #kkt_solver.kkt_err_norm.ratio > 1e-1
        @show kkt_solver.kkt_err_norm
        return :failure, iter, Class_ls_info()
    end

    success, iter, step_info = simple_ls(iter, kkt_solver.dir, ls_mode, filter, pars, min_step_size, timer)

    return success, iter, step_info
end
