
type Class_homog <: abstract_nlp
    nlp::abstract_nlp
    fshift::Float64

    # min f(x) + lambda_vec * (x - centre_point)^2

    function Class_homog(nlp::abstract_nlp)
        return new(nlp)
    end
end

function convert_to_homog!(iter::Class_iterate, pars::Class_parameters)

    iter.nlp = Class_homog(iter.nlp)
    iter.nlp.fshift = 1000.0
    p = iter.point
    p.x = [p.x; 1.0];
    p.s = [p.s; 1.0];
    p.y = [p.y; p.mu];
    iter.primal_residual_intial = [iter.primal_residual_intial; 0.0]
end

function x_true(x::Array{Float64,1})
    return x[1:(end - 1)] / x[end]
end

function eval_f(m::Class_homog, x::Array{Float64,1})
    return x[end] * (eval_f(m.nlp, x_true(x)) - m.fshift)
end

function eval_a(m::Class_homog, x::Array{Float64,1})
    tau = x[end]
    return [tau * eval_a(m.nlp, x_true(x)); tau];
end


function eval_jac(m::Class_homog, x::Array{Float64,1})
    #@show "eval_jac"
    xhat = x_true(x)
    J = eval_jac(m.nlp, xhat)
    b = eval_a(m.nlp, xhat) - J * xhat

    return [[J b]; [zeros(length(xhat))' 1.0]];
end

function eval_grad_lag(m::Class_homog, x::Array{Float64,1}, y::Array{Float64,1}, w::Float64=1.0)
    #@show "eval_grad_lag"
    xhat = x_true(x)
    g = eval_grad_lag(m.nlp, xhat, zeros(length(y)-1), w)
    #a = eval_a(m.nlp, xhat)
    df_dtau = (eval_f(m.nlp, xhat) - m.fshift) - dot(xhat, g)

    return [g; df_dtau] - eval_jac(m, x)' * y;
end


function eval_lag_hess(m::Class_homog, x::Array{Float64,1}, y::Array{Float64,1}, w::Float64)
    #@show "eval_grad_hess"
    xhat = x_true(x)
    H = eval_lag_hess(m.nlp, xhat,  y[1:(end-1)], w);
    M = [[H  H * xhat];
    [xhat' * H  dot(xhat, H * xhat)]];

    return  M / y[end]
end
