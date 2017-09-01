
type Class_proximal <: abstract_nlp
    nlp::abstract_nlp
    #lambda_vec::Array{Float64,1}
    #centre_point::Array{Float64,1}
    lambda::Float64
    style::Symbol

    x_norm_penalty::Float64
    a_norm_penalty::Float64

    # min f(x) + lambda_vec * (x - centre_point)^2
    function Class_proximal(nlp::abstract_nlp, pars::Class_parameters)
        return new(nlp, 0.0, :sqrt, pars.x_norm_penalty, pars.a_norm_penalty)
    end

    #function Class_proximal(nlp::abstract_nlp, lambda_vec::Array{Float64,1}, centre_point::Array{Float64,1})
    #    return new(nlp, lambda_vec, centre_point, :sqrt)
    #end

    #function Class_proximal(nlp::abstract_nlp, lambda_vec::Array{Float64,1}, centre_point::Array{Float64,1})
    #    return new(nlp, lambda_vec, centre_point, :sqrt)
    #end
    #
    #function Class_proximal(nlp::abstract_nlp, lambda_vec::Array{Float64,1})
    #    return Class_proximal(nlp, lambda_vec, zeros(length(lambda_vec)))
    #end
end

function convert_to_prox!(iter::Class_iterate, pars::Class_parameters, lambda::Float64)
    iter.nlp = Class_proximal(iter.nlp, pars)
    update_prox!(iter, pars, lambda)
end

function update_prox!(iter::Class_iterate, pars::Class_parameters, lambda::Float64)
    prox = iter.nlp::Class_proximal
    prox.lambda = lambda
end

function eval_f(m::Class_proximal, x::Array{Float64,1})
    if m.style == :quad
      error("")
      #return eval_f(m.nlp, x) + 0.5 * sum(m.lambda_vec .* (x - m.centre_point).^2)
    elseif m.style == :sqrt
      β = m.x_norm_penalty
      x_norm_pen = β * sum(sqrt(x.^2 + 1.0 / β^2))
      a_norm_pen = m.a_norm_penalty * sum(eval_a(m.nlp,x))
      rx = x_norm_pen + a_norm_pen
      return eval_f(m.nlp, x) + m.lambda * rx
    end
end

function eval_a(m::Class_proximal, x::Array{Float64,1})
    return eval_a(m.nlp, x)
end

function eval_jac(m::Class_proximal, x::Array{Float64,1})
    return eval_jac(m.nlp, x)
end

function eval_grad_f(m::Class_proximal, x::Array{Float64,1})
    if m.style == :quad
        return eval_grad_f(m.nlp, x) + m.lambda * x
    elseif m.style == :sqrt
        β = m.x_norm_penalty
        x_norm_grad =  β * x ./ sqrt(x.^2 + 1.0 / β^2)
        a_norm_grad = m.a_norm_penalty * eval_jac(m.nlp,x)' * ones(length(eval_a(m.nlp,x)))
        rx_g = x_norm_grad + a_norm_grad
        return eval_grad_f(m.nlp, x) + m.lambda * rx_g
    end
end

function eval_Jt_prod(m::Class_proximal, x::Array{Float64,1}, y::Array{Float64,1})
    return eval_Jt_prod(m.nlp, x, y)
end

function eval_lag_hess(m::Class_proximal, x::Array{Float64,1}, y::Array{Float64,1}, w::Float64)
    if m.style == :quad
      return eval_lag_hess(m.nlp, x, y, w) + spdiagm(m.lambda_vec)
    elseif m.style == :sqrt
      y_rx = m.lambda * m.a_norm_penalty

      β = m.x_norm_penalty
      denominator = (β^2 * x.^2 + 1.0) .* sqrt( x.^2 + 1.0 / β^2 )
      H_x_norm = spdiagm(m.lambda ./ denominator)
      return eval_lag_hess(m.nlp, x, y + y_rx, w) + H_x_norm
    end
    #return eval_lag_hess(m.nlp, x, y, w)
end
