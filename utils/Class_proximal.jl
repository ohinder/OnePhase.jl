
type Class_proximal <: abstract_nlp
    nlp::abstract_nlp
    lambda_vec::Array{Float64,1}
    centre_point::Array{Float64,1}

    # min f(x) + lambda_vec * (x - centre_point)^2

    function Class_proximal(nlp::abstract_nlp, lambda_vec::Array{Float64,1}, centre_point::Array{Float64,1})
        return new(nlp, lambda_vec, centre_point)
    end

    function Class_proximal(nlp::abstract_nlp, lambda_vec::Array{Float64,1})
        return Class_proximal(nlp, lambda_vec, zeros(length(lambda_vec)))
    end
end

function convert_to_prox!(iter::Class_iterate, pars::Class_parameters)
    iter.nlp = Class_proximal(iter.nlp, zeros(length(get_x(iter))))
    update_prox!(iter, pars)
end

function update_prox!(iter::Class_iterate, pars::Class_parameters)
    prox = iter.nlp::Class_proximal
    x = get_x(iter)
    #prox.lambda_vec = min(100.0, get_mu(iter)) * ones(length(x)) / (10.0 * (1.0 + norm(x, Inf) ) )
    #prox.lambda_vec = get_mu(iter) * ones(length(x)) ./ (norm(x, 2) * 10.0 + 100.0)
    if pars.proximal_style == :none
      lambda = 0.0
    elseif pars.proximal_style == :test1
      lambda = min(100.0, get_mu(iter)) / ((norm(x, Inf)^2 + 1.0) * 10.0)
    elseif pars.proximal_style == :test2
      lambda = min(100.0, get_mu(iter)) / ((norm(x, Inf) + 1.0) * 10.0)
    elseif pars.proximal_style == :test3
      lambda = min(1.0, get_mu(iter)) / ((norm(x, Inf) + 1.0) * 10.0)
    end
    prox.lambda_vec = lambda * ones(length(x))
end


#=
function nbounds_orginal(nlp::Class_CUTEst)
    return length(nlp.bvar.l_i) + length(nlp.bvar.u_i)
end

function ncons_orginal(nlp::Class_CUTEst)
    return length(nlp.bcon.l_i) + length(nlp.bcon.u_i)
end

function cons_indicies(nlp::Class_CUTEst)
    m = ncons_orginal(nlp)
    if m > 0
      return 1:m
    else
      return []
    end
end


function bound_indicies(nlp::Class_CUTEst)
    m = ncons_orginal(nlp)
    r = nbounds_orginal(nlp)
    if r > 0
      return (m + 1):(m + r)
    else
      return [];
    end
end

function y_l_con(y::Array{Float64,1}, m::Class_CUTEst)
    return y[1:length(m.bcon.l)]
end

function y_u_con(y::Array{Float64,1}, m::Class_CUTEst)
    n_lcon = length(m.bcon.l)
    return y[(n_lcon + 1):(n_lcon + length(m.bcon.u))]
end=#

function eval_f(m::Class_proximal, x::Array{Float64,1})
    return eval_f(m.nlp, x) + 0.5 * sum(m.lambda_vec .* (x - m.centre_point).^2)
end

function eval_a(m::Class_proximal, x::Array{Float64,1})
    return eval_a(m.nlp, x)
end


function eval_jac(m::Class_proximal, x::Array{Float64,1})
    return eval_jac(m.nlp, x)
end

function eval_grad_f(m::Class_proximal, x::Array{Float64,1})
    return eval_grad_f(m.nlp, x) + m.lambda_vec .* (x - m.centre_point)
end

function eval_Jt_prod(m::Class_proximal, x::Array{Float64,1}, y::Array{Float64,1})
    return eval_Jt_prod(m.nlp, x, y)
end

function eval_lag_hess(m::Class_proximal, x::Array{Float64,1}, y::Array{Float64,1}, w::Float64)
    return eval_lag_hess(m.nlp, x, y, w) + spdiagm(m.lambda_vec)
end
