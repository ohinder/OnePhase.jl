
type Class_QP <: abstract_nlp
    b::Array{Float64,1}
    A::SparseMatrixCSC{Float64,Int64}
    c::Array{Float64,1}
    Q::SparseMatrixCSC{Float64,Int64}
    # A x >= b
    function Class_QP(b::Array{Float64,1}, A::SparseMatrixCSC{Float64,Int64}, c::Array{Float64,1}, Q::SparseMatrixCSC{Float64,Int64})
        return new(b, A, c, Q)
    end

    function Class_QP(b::Array{Float64,1}, A::SparseMatrixCSC{Float64,Int64}, c::Array{Float64,1})
        n = length(c)
        return new(b, A, c, spzeros(n,n))
    end
end

function dim(QP::Class_QP)
  return length(QP.c)
end

function ncon(QP::Class_QP)
  return length(QP.b)
end

function eval_f(QP::Class_QP, x::Array{Float64,1})
    return dot(QP.c, x) + 0.5 * dot(x, QP.Q * x)
end

function eval_a(QP::Class_QP, x::Array{Float64,1})
    return QP.A * x - QP.b
end

function eval_jac(QP::Class_QP, x::Array{Float64,1})
    return QP.A
end

function eval_grad_lag(QP::Class_QP, x::Array{Float64,1}, y::Array{Float64,1}, w::Float64)
    return w * (QP.c + QP.Q * x) - QP.A' * y
end

function eval_lag_hess(QP::Class_QP, x::Array{Float64,1}, y::Array{Float64,1}, w::Float64)
    return w * QP.Q
end

function Construct_LP_problem(A::SparseMatrixCSC{Float64,Int64}, b::Array{Float64,1}, c::Array{Float64,1})
  # min c'x
  # s.t. A x >= b
  LP = Class_nlp()
  sigma = 0.0#1e-4

  LP.dim = length(c)
  LP.ncon = length(b)
  LP.eval_f = x -> dot(c, x) + sigma * dot(x, x);
  LP.eval_a = x -> (A * x - b);
  LP.eval_jac = x -> A;
  LP.eval_grad_lag = (x, y, w) ->  w * c + 2 * sigma * x - A' * y;
  LP.eval_lag_hess = (x, y, w) -> 2 * sigma * speye(length(c),length(c));

  return LP
end
