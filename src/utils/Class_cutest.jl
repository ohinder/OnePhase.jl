importall CUTEst
#using SparseArrays

export get_original_x, get_y

mutable struct Class_bounds
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



function _i_not_fixed(m::NLPModels.AbstractNLPModel)
    return (1:m.meta.nvar)[m.meta.lvar .!= m.meta.uvar]
end

mutable struct Class_CUTEst <: abstract_nlp
    nlp::NLPModels.AbstractNLPModel

    bcon::Class_bounds
    bvar::Class_bounds

    function Class_CUTEst(nlp::NLPModels.AbstractNLPModel)
        ind = _i_not_fixed(nlp)
        return new(nlp, Class_bounds(nlp.meta.lcon[:], nlp.meta.ucon[:]), Class_bounds(nlp.meta.lvar[ind], nlp.meta.uvar[ind]))
    end
end

#=
function con_info(m::Class_CUTEst, index::Int64)
    con = false
    if index in cons_indicies(m)
      con = true
      if index < length(nlp.bcon.l_i)
        #[lb(a, m.bcon); ub(a, m.bcon); lb(x, m.bvar); ub(x, m.bvar)]
      end
    end



    return;
end=#

function suggested_starting_point(m::Class_CUTEst)
    ind = _i_not_fixed(m.nlp)
    return deepcopy(m.nlp.meta.x0[ind])
end

function ncon(m::Class_CUTEst)
    return nbounds_orginal(m) + ncons_orginal(m)
end

function lb(x::Array{Float64,1}, bd::Class_bounds)
    return x[bd.l_i] - bd.l
end

function ub(x::Array{Float64,1}, bd::Class_bounds)
    return bd.u - x[bd.u_i]
end

function linear_cons(m::Class_CUTEst)
    is_linear = zeros(m.nlp.meta.ncon)
    is_linear[m.nlp.meta.lin] = 1.0
    vec = [is_linear[m.bcon.l_i]; is_linear[m.bcon.u_i]; ones(nbounds_orginal(m))]

    return 1.0 .== vec
end

function ineq_cons(m::Class_CUTEst)
    is_ineq = zeros(m.nlp.meta.ncon)
    is_ineq[m.nlp.meta.lcon .== m.nlp.meta.ucon] = 1.0
    vec = [is_ineq[m.bcon.l_i]; is_ineq[m.bcon.u_i]; ones(nbounds_orginal(m))]

    return 1.0 .== vec
end

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

#function linear_indicies(nlp::Class_CUTEst)
#    #nlp.nlp.meta.lin
#end


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
end

function eval_f(m::Class_CUTEst, x::Array{Float64,1})
    return obj(m.nlp, _cute_x(m, x) )
end

function eval_a(m::Class_CUTEst, x::Array{Float64,1})
    a = cons(m.nlp, _cute_x(m, x) )
    return [lb(a, m.bcon); ub(a, m.bcon); lb(x, m.bvar); ub(x, m.bvar)];
end

function _cute_x(m::Class_CUTEst, x::Array{Float64,1})
    ind = _i_not_fixed(m.nlp)
    #@show length(ind)
    if(length(x) != length(ind))
      error("$(length(x)) = length(x) != length(i_not_fixed) = $(length(ind))")
    end


    cute_x = deepcopy(m.nlp.meta.lvar) # get correct values of fixed variables
    cute_x[ind] = x

    return cute_x
end

function eval_jac(m::Class_CUTEst, x::Array{Float64,1})
    cute_x = _cute_x(m, x)
    J_full_T = jac(m.nlp, cute_x)'
    J_T = J_full_T[_i_not_fixed(m.nlp),:];
    my_eye = SparseArrays.speye(length(x))
    Q_T = [J_T[:,m.bcon.l_i] -J_T[:,m.bcon.u_i] my_eye[:,m.bvar.l_i] -my_eye[:,m.bvar.u_i]];
    return Q_T'

    #return @time [J[m.bcon.l_i,:]; -J[m.bcon.u_i,:]; my_eye[m.bvar.l_i,:]; -my_eye[m.bvar.u_i,:]];
end

function eval_grad_f(m::Class_CUTEst, x::Array{Float64,1})
    return grad(m.nlp, _cute_x(m, x))[_i_not_fixed(m.nlp)]
end

function get_constrduals(m::Class_CUTEst, y::Array{Float64,1})
     constrdual = zeros(m.nlp.meta.ncon)
     st_u = length(m.bcon.l) #+ ncons_orginal(m)
     constrdual[m.bcon.l_i] += y[1:st_u]
     constrdual[m.bcon.u_i] -= y[(st_u+1):(length(m.bcon.u) + st_u)]
     return constrdual
end

function get_reducedcosts(m::Class_CUTEst, y::Array{Float64,1})
    rc = zeros(m.nlp.meta.nvar)
    st_l = ncons_orginal(m)
    st_u = length(m.bvar.l) + ncons_orginal(m)
    rc[m.bvar.l_i] += y[(st_l+1):st_u]
    rc[m.bvar.u_i] -= y[(st_u+1):(length(m.bvar.u) + st_u)]
    return rc
end

function y_cons_net(m::Class_CUTEst, y::Array{Float64,1})
    return;
end

function eval_jtprod(m::Class_CUTEst, x::Array{Float64,1}, y::Array{Float64,1})
    return jtprod(nlp_raw, x, y)
end

function eval_schur_and_jtprod(m::Class_CUTEst, x::Array{Float64,1}, y::Array{Float64,1}, s::Array{Float64,1})

end

#=function ncon(m::Class_CUTEst)
    return length(m.bcon.l) + length(m.bcon.u)
end=#

#=function nvar(m::Class_CUTEst)
    return m.nlp.meta.nvar
end=#
#=
function make_symmetric(M::SparseMatrixCSC{Float64,Int32})
    n = size(M,1)
    new_M = spzeros(n,n)
    rows = rowvals(M)
    vals = nonzeros(M)
    for col = 1:n
      for j in nzrange(M, col)
         row = rows[j]
         val = vals[j]
         # perform sparse wizardry...
         new_M[col,row] = vals[j]
         new_M[row,col] = vals[j]
      end
    end

    return new_M
end=#

function eval_lag_hess(m::Class_CUTEst, x::Array{Float64,1}, y::Array{Float64,1}, w::Float64)
    y_cons = zeros(m.nlp.meta.ncon)
    y_cons[m.bcon.l_i] -= y_l_con(y, m)
    y_cons[m.bcon.u_i] += y_u_con(y, m)

    H = hess(m.nlp, _cute_x(m, x), obj_weight=w, y=y_cons);
    #H = hess(m.nlp, _cute_x(m, x), w, y_cons);
    #@show typeof(H)

    ind = _i_not_fixed(m.nlp)

    H_not_fixed = H[ind,ind]
    #H_true = SparseArrays.sparse(Symmetric(H_not_fixed,:L) + spzeros(length(x),length(x)))
    #H_true = H_not_fixed + H_not_fixed' - spdiagm(diag(H_not_fixed))
    #convert()
    #@time H_true = make_symmetric(H_not_fixed)
    #@show LinearAlgebra.norm(H_true - H2,Inf)

    return H_not_fixed
end

function eval_Jt_prod(m::Class_CUTEst, x::Array{Float64,1}, y::Array{Float64,1})
    return jtprod(m.nlp, x, y)
end
