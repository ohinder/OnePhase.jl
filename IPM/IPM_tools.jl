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
    return norm(eval_grad_lag(iter, iter.point.mu),Inf) * dual_scale(iter, pars)
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
