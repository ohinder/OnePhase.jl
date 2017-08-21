
type Class_proximal <: abstract_nlp
    nlp::abstract_nlp
    lambda_vec::Array{Float64,1}
    centre_point::Array{Float64,1}
    style::Symbol
    beta::Float64

    # min f(x) + lambda_vec * (x - centre_point)^2

    function Class_proximal(nlp::abstract_nlp, lambda_vec::Array{Float64,1}, centre_point::Array{Float64,1})
        return new(nlp, lambda_vec, centre_point, :sqrt, 1e-4)
    end

    function Class_proximal(nlp::abstract_nlp, lambda_vec::Array{Float64,1})
        return Class_proximal(nlp, lambda_vec, zeros(length(lambda_vec)))
    end
end

function convert_to_prox!(iter::Class_iterate, pars::Class_parameters, lambda::Float64)
    iter.nlp = Class_proximal(iter.nlp, zeros(length(get_x(iter))))
    update_prox!(iter, pars, lambda)
end

function update_prox!(iter::Class_iterate, pars::Class_parameters, lambda::Float64)
    prox = iter.nlp::Class_proximal
    x = get_x(iter)
    #prox.lambda_vec = min(100.0, get_mu(iter)) * ones(length(x)) / (10.0 * (1.0 + norm(x, Inf) ) )
    #prox.lambda_vec = get_mu(iter) * ones(length(x)) ./ (norm(x, 2) * 10.0 + 100.0)
    #=if pars.proximal_style == :none
      lambda = 0.0
    elseif pars.proximal_style == :test1
      lambda = min(100.0, get_mu(iter)) / ((norm(x, Inf)^2 + 1.0) * 10.0)
    elseif pars.proximal_style == :test2
      lambda = min(100.0, get_mu(iter)) / ((norm(x, Inf) + 1.0) * 10.0)
    elseif pars.proximal_style == :test3
      lambda = min(1.0, get_mu(iter)) / ((norm(x, Inf) + 1.0) * 10.0)
    elseif pars.proximal_style == :fixed
      lambda = get_mu(iter)
    end=#
    prox.lambda_vec = lambda * ones(length(x))
end

beta_thing = -1e-4 #-1e-8

function eval_f(m::Class_proximal, x::Array{Float64,1})
    if m.style == :quad
      return eval_f(m.nlp, x) + 0.5 * sum(m.lambda_vec .* (x - m.centre_point).^2)
    elseif m.style == :sqrt
      return eval_f(m.nlp, x) + sum(m.beta * m.lambda_vec .* sqrt((x - m.centre_point).^2 + 1.0 / m.beta^2)) - norm(m.lambda_vec,Inf) * beta_thing * sum(eval_a(m.nlp,x))
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
        return eval_grad_f(m.nlp, x) + m.lambda_vec .* (x - m.centre_point)
    elseif m.style == :sqrt
      return eval_grad_f(m.nlp, x) + m.lambda_vec * m.beta .* x ./ sqrt((x - m.centre_point).^2 + 1.0 / m.beta^2) - norm(m.lambda_vec,Inf) * beta_thing * eval_jac(m.nlp,x)' * ones(length(eval_a(m.nlp,x)))
    end
end

function eval_Jt_prod(m::Class_proximal, x::Array{Float64,1}, y::Array{Float64,1})
    return eval_Jt_prod(m.nlp, x, y)
end

function eval_lag_hess(m::Class_proximal, x::Array{Float64,1}, y::Array{Float64,1}, w::Float64)
    if m.style == :quad
      return eval_lag_hess(m.nlp, x, y, w) + spdiagm(m.lambda_vec)
    elseif m.style == :sqrt
      denominator = (m.beta^2 * x.^2 + 1.0) .* sqrt((x - m.centre_point).^2 + 1.0 / m.beta^2)
      return eval_lag_hess(m.nlp, x, y -  norm(m.lambda_vec,Inf) * beta_thing, w) + spdiagm(m.beta * m.lambda_vec ./ denominator)
    end
    #return eval_lag_hess(m.nlp, x, y, w)
end
