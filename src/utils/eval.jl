using SparseArrays

function comp_merit(it::Class_iterate)
    return LinearAlgebra.norm(comp(it),Inf)^3 / get_mu(it)^2
end

#function comp_merit_predicted(it::Class_iterate, dir::Class_point, step_size::Float64)
#    return comp_predicted(it, dir, step_size)^3
#end

function comp(it::Class_iterate)
    return it.point.s .* it.point.y .- it.point.mu
end



function comp_ratio_max(it::Class_iterate)
    return max(maximum(it.point.s .* it.point.y ./ it.point.mu), maximum(it.point.mu ./ (it.point.s .* it.point.y)))
end

mutable struct Eval_NaN_error <: Exception
   num::Float64
   point::Array{Float64,1}
   place::String

   function Eval_NaN_error(num::Float64,point::Array{Float64,1},place::String)
     return new(num::Float64,point::Array{Float64,1},place::String)
   end
end

Base.showerror(io::IO, e::Eval_NaN_error) = print(io, "$(e.num) by $(e.place) not expected!");

function isbad(num::Float64)
    if isnan(num) || isinf(num)
      return true
      #my_warn("$num not expected!")
      #Cutest_NaN_error("$num not expected!", num, point, place)
    else
      return false
    end
end

function isbad(vec)
  for i = 1:length(vec)
    if isbad(vec[i])
      return true
    end
  end

  return false
end

function eval_r(it::Class_iterate)
    x = it.point.x
    a_norm_pen = it.a_norm_penalty_par * sum(get_cons(it))
    return a_norm_pen
end

function eval_grad_r(it::Class_iterate)
    x = it.point.x
    a_norm_grad = it.a_norm_penalty_par * eval_jac_T_prod(it, ones(length(get_cons(it))))
    return a_norm_grad
end

function get_fval(it::Class_iterate)
    return it.cache.fval
end

function get_grad(it::Class_iterate)
    return it.cache.grad #+
end

function get_cons(it::Class_iterate)
    return it.cache.cons
end

function get_jac(it::Class_iterate)
    return it.cache.J
end

function get_jac_T(it::Class_iterate)
    return it.cache.J_T
end

function eval_J_T_J(it::Class_iterate, diag_vals::Vector)
    #return it.cache.J_T * spdiagm(diag_vals) * it.cache.J
    return it.cache.J_T * sparse(Diagonal(diag_vals)) * it.cache.J
end

function eval_diag_J_T_J(it::Class_iterate,diag_vals::Vector)
    n = dim(it)
    di = zeros(n)
    for i = 1:n
        a = it.cache.J[:,i]
        for j in a.nzind
            di[i] += a[j]^2 * diag_vals[j]
        end
        #dot(a .* diag_vals,a)
    end
    return di
end

function eval_jac_prod(it::Class_iterate, x::Vector)
    return it.cache.J * x
end

function eval_jac_T_prod(it::Class_iterate, y::Vector)
    return it.cache.J_T * y
end

function get_primal_res(it::Class_iterate)
    return it.cache.cons .- it.point.s
end

function get_max_vio(it::Class_iterate)
    return -min(0.0,minimum(it.cache.cons))
end

function eval_phi(it::Class_iterate, mu::Float64)
    if all(it.point.s .> 0.0)
        return get_fval(it) - mu * sum( log.( it.point.s ) ) + mu * eval_r(it)
    else
        return Inf
    end
end

function eval_grad_phi(it::Class_iterate, mu::Float64)
    y_tilde = it.point.mu ./ it.point.s #
    return eval_grad_lag(it, mu, y_tilde)
end


function get_lag_hess(it::Class_iterate)
    return it.cache.H # remember this is a triangular matrix!!!!
end

function eval_grad_lag(it::Class_iterate, mu::Float64)
    return eval_grad_lag(it, mu, it.point.y)
end

function eval_grad_lag(it::Class_iterate, mu::Float64, y::Array{Float64,1})
    return get_grad(it) - eval_jac_T_prod(it, y) + mu * eval_grad_r(it)
end

function dynamic_eval_Jt_prod(it::Class_iterate)
    return eval_Jt_prod(it.nlp, it.point.x, it.point.y)
end


function eval_farkas(it::Class_iterate, y::Array{Float64,1})
    return max(farkas_certificate(it, y), stationary_infeasible_measure(it,y))
end

function farkas_certificate(it::Class_iterate, y::Array{Float64,1})
    c = it.cache.cons
    feas_obj = -sum(c .* y)
    if feas_obj > 0.0
      return LinearAlgebra.norm(it.cache.J_T * y, 1) / feas_obj
    else
      return Inf
    end
end


function stationary_infeasible_measure(it::Class_iterate, y::Array{Float64,1})
    return (LinearAlgebra.norm(it.cache.J_T * y, 1) + dot(it.point.s, it.point.y)) / LinearAlgebra.norm(y,1)
end

function eval_farkas(it::Class_iterate)
    return eval_farkas(it, it.point.y)
end

function eval_merit_function(it::Class_iterate, pars::Class_parameters)
    if is_feasible(it, pars.ls.comp_feas)
        if length(it.point.s) > 0
          comp_penalty = LinearAlgebra.norm(comp(it), Inf)^3 / (it.point.mu)^2
        else
          comp_penalty = 0.0
        end

        return eval_phi(it, it.point.mu) + comp_penalty
    else
        return Inf
    end
end

function eval_phi_diff(it::Class_iterate, candidate::Class_iterate, mu::Float64)
    fdiff = get_fval(candidate) - get_fval(it)
    rdiff = mu * eval_r(candidate) - mu * eval_r(it)
    logdiff = - mu * sum( log.( candidate.point.s ./ it.point.s) )

    return fdiff + rdiff + logdiff
end

function eval_merit_function_difference(it::Class_iterate, candidate::Class_iterate, pars::Class_parameters)
    if is_feasible(candidate, pars.ls.comp_feas)
        if length(it.point.s) > 0
          comp_penalty = LinearAlgebra.norm(comp(candidate), Inf)^3 / (it.point.mu)^2 - LinearAlgebra.norm(comp(it), Inf)^3 / (it.point.mu)^2
        else
          comp_penalty = 0.0
        end

        #phi_diff = eval_phi(candidate, candidate.point.mu) - eval_phi(it, it.point.mu)
        phi_diff = eval_phi_diff(it, candidate, candidate.point.mu)
        return phi_diff + comp_penalty
    else
        return Inf
    end
end

function is_lower_triangular(mat)
    for i = 1:size(mat,1)
        for j = (i+1):size(mat,2)
            if mat[i,j] != 0.0
                return false
            end
        end
    end
    return true
end

function vector_product(LowerTriangularMat::SparseMatrixCSC{Float64,Int64}, vector::Array{Float64,1})
    DEBUG_MODE = false
    if DEBUG_MODE
        warn("checking lower triangular property")
        @assert is_lower_triangular(LowerTriangularMat) == true
    end
    v1 = LowerTriangularMat * vector
    v2 = LowerTriangularMat' * vector
    return v1 + v2 - diag(LowerTriangularMat) .* vector
end

function hess_product(it::Class_iterate, vector::Array{Float64,1})
    return vector_product(it.cache.H, vector)
end

function phi_predicted_reduction_primal_dual(it::Class_iterate, dir::Class_point, step_size::Float64)
    H = it.cache.H
    J = it.cache.J

    v = J * dir.x
    J_gain = dot(v.^2, it.point.y ./ it.point.s)

    #@show H

    h_prod = vector_product(H, dir.x) # the hessian is lower triangular!

    return step_size * dot(dir.x, eval_grad_phi(it, it.point.mu)) + step_size^2 * 0.5 * (dot(dir.x, h_prod) + J_gain)
end

function comp_predicted(it::Class_iterate, dir::Class_point, step_size::Float64)
    y = it.point.y;
    s = it.point.s;
    return s .* y + dir.y .* s * step_size + dir.s .* y * step_size .- (it.point.mu + dir.mu * step_size)
end

function merit_function_predicted_reduction(it::Class_iterate, dir::Class_point, step_size::Float64)
    C_k = LinearAlgebra.norm(comp(it),Inf)
    P_k = LinearAlgebra.norm(comp_predicted(it, dir, step_size), Inf)
    #@show C_k, P_k
    if length(it.point.s) > 0
      comp_penalty = (P_k^3 - C_k^3) / (it.point.mu)^2
    else
      comp_penalty = 0.0
    end
    #@show comp_penalty

    phi_red = phi_predicted_reduction_primal_dual(it, dir, step_size)
    #@show predict_red

    #@show comp_penalty

    return phi_red + comp_penalty
end

function eval_kkt_err(it::Class_iterate, par::Class_parameters)
    return scaled_dual_feas(it, par) + LinearAlgebra.norm(comp(it), Inf)
end
