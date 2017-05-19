using CUTEst

type Class_bounds
    l_i::Array{Int64,1}
    u_i::Array{Int64,1}
    l::Array{Float64,1}
    u::Array{Float64,1}

    function Class_bounds(lb_vec::Array{Float64,1}, ub_vec::Array{Float64,1})
        @assert(length(lb_vec) == length(ub_vec))
        this = new([],[],[],[]);
        for i = 1:length(lb_vec)
            if lb_vec[i] > -Inf
               push!(this.l_i, i)
               push!(this.l, lb_vec[i])
            end

            if ub_vec[i] < Inf
               push!(this.u_i, i)
               push!(this.u, ub_vec[i])
            end
        end

        return this;
    end
end


type Class_CUTEst <: abstract_nlp
    nlp::CUTEst.CUTEstModel

    bcon::Class_bounds
    bvar::Class_bounds

    function Class_CUTEst(nlp::CUTEst.CUTEstModel)
        return new(nlp, Class_bounds(nlp.meta.lcon, nlp.meta.ucon), Class_bounds(nlp.meta.lvar, nlp.meta.uvar))
    end
end


function lb(x::Array{Float64,1}, bd::Class_bounds)
    return x[bd.l_i] - bd.l
end

function ub(x::Array{Float64,1}, bd::Class_bounds)
    return bd.u - x[bd.u_i]
end

function y_l_con(y::Array{Float64,1}, m::Class_CUTEst)
    return y[1:length(m.bcon.l)]
end

function y_u_con(y::Array{Float64,1}, m::Class_CUTEst)
    n_lcon = length(m.bcon.l)
    return y[(n_lcon + 1):(n_lcon + length(m.bcon.u))]
end

function eval_f(m::Class_CUTEst, x::Array{Float64,1})
    return obj(m.nlp, x)
end

function eval_a(m::Class_CUTEst, x::Array{Float64,1})
    a = cons(m.nlp, x)
    return [lb(a, m.bcon); ub(a, m.bcon); lb(x, m.bvar); ub(x, m.bvar)];
end

function eval_jac(m::Class_CUTEst, x::Array{Float64,1})
    J = jac(m.nlp, x);
    my_eye = speye(length(x))
    return [J[m.bcon.l_i,:]; -J[m.bcon.u_i,:]; my_eye[m.bvar.l_i,:]; -my_eye[m.bvar.u_i,:]];
end

function eval_grad_lag(m::Class_CUTEst, x::Array{Float64,1}, y::Array{Float64,1}, w::Float64)
    #y_cons = m.nlp.meta.lcon;
    #y_vars =
    J = eval_jac(m, x)
    #@show size(J), length(x), length(y)
    return w * grad(m.nlp, x) - J' * y
    #+ jprod(m, x, y_cons) + y_vars
end

function ncon(m::Class_CUTEst)
    return m.nlp.meta.ncon
end

function nvar(m::Class_CUTEst)
    return m.nlp.meta.nvar
end


function eval_lag_hess(m::Class_CUTEst, x::Array{Float64,1}, y::Array{Float64,1}, w::Float64)
    y_cons = zeros(ncon(m))
    y_cons[m.bcon.l_i] -= y_l_con(y, m)
    y_cons[m.bcon.u_i] += y_u_con(y, m)
    return hess(m.nlp, x, obj_weight=w, y=y_cons)
end
