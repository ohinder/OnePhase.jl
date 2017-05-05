using CUTEst


type Class_CUTEst <: abstract_nlp
    nlp::CUTEst.CUTEstModel

    function Class_CUTEst(nlp::CUTEst.CUTEstModel)
        return new(nlp)
    end
end

function dim(m::Class_CUTEst)
  return m.meta.nvar
end

function ncon(m::Class_CUTEst)
  return m.meta.ncon + m.meta.nvar
end

function eval_f(m::Class_CUTEst, x::Array{Float64,1})
    return obj(m, x)
end

function eval_a(m::Class_CUTEst, x::Array{Float64,1})
    return cons(m, x)
end

function eval_jac(m::Class_CUTEst, x::Array{Float64,1})
    return jac(m, x)
end

function eval_grad_lag(m::Class_CUTEst, x::Array{Float64,1}, y::Array{Float64,1}, w::Float64)
    total_cons = m.meta.ncon

    y_true[m.meta.la]
    return w * grad(m, x) + jprod(m, x, y_true)
end

function eval_lag_hess(m::Class_CUTEst, x::Array{Float64,1}, y::Array{Float64,1}, w::Float64)
    return w * CUTEst.Q
end
