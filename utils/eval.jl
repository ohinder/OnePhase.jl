function comp_merit(it::Class_iterate)
    return norm(comp(it),Inf)^3 / get_mu(it)^2
end

#function comp_merit_predicted(it::Class_iterate, dir::Class_point, step_size::Float64)
#    return comp_predicted(it, dir, step_size)^3
#end

function comp(it::Class_iterate)
    return it.point.s .* it.point.y - it.point.mu
end



function comp_ratio_max(it::Class_iterate)
    return max(maximum(it.point.s .* it.point.y ./ it.point.mu), maximum(it.point.mu ./ (it.point.s .* it.point.y)))
end

type Eval_NaN_error <: Exception
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
    beta = it.x_norm_penalty_par
    x_norm_pen = beta * sum(sqrt(x.^2 + 1.0 / beta^2))
    a_norm_pen = it.a_norm_penalty_par * sum(get_cons(it))
    return (x_norm_pen + a_norm_pen) * it.use_reg
end

function eval_grad_r(it::Class_iterate)
    x = it.point.x
    beta = it.x_norm_penalty_par
    #@show it.x_norm_penalty_par
    x_norm_grad =  beta * x ./ sqrt(x.^2 + 1.0 / beta^2)
    a_norm_grad = it.a_norm_penalty_par * eval_jac_T_prod(it, ones(length(get_cons(it))))
    return (x_norm_grad + a_norm_grad) * it.use_reg
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

#function get_jacobian(it::Class_iterate)
#    return it.cache.J # should make this smaller and use the sparse non-zero format
#end

function eval_J_T_J(it::Class_iterate, diag_vals::Vector)
    return it.cache.J_T * spdiagm(diag_vals) * it.cache.J
end

function eval_jac_prod(it::Class_iterate, x::Vector)
    return it.cache.J * x
end

function eval_jac_T_prod(it::Class_iterate, y::Vector)
    return it.cache.J_T * y
end

function get_primal_res(it::Class_iterate)
    return it.cache.cons - it.point.s
end

function get_max_vio(it::Class_iterate)
    return -min(0.0,minimum(it.cache.cons))
end

function eval_phi(it::Class_iterate, mu::Float64)
    return get_fval(it) - mu * sum( log( it.point.s ) ) + mu * eval_r(it)
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
    feas_obj = -mean(it.cache.cons .* y)
    if feas_obj > 0.0
      return norm(y' * it.cache.J, Inf) /  min(norm(y, Inf), feas_obj)
    else
      return Inf
    end
end

function eval_farkas(it::Class_iterate)
    return eval_farkas(it, it.point.y)
end

function eval_merit_function(it::Class_iterate, pars::Class_parameters)
    if is_feasible(it, pars.comp_feas)
        if length(it.point.s) > 0
          comp_penalty = norm(comp(it), Inf)^3 / (it.point.mu)^2
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
    logdiff = - mu * sum( log( candidate.point.s ./ it.point.s) )

    return fdiff + rdiff + logdiff
end

function eval_merit_function_difference(it::Class_iterate, candidate::Class_iterate, pars::Class_parameters)
    if is_feasible(candidate, pars.comp_feas)
        if length(it.point.s) > 0
          comp_penalty = norm(comp(candidate), Inf)^3 / (it.point.mu)^2 - norm(comp(it), Inf)^3 / (it.point.mu)^2
        else
          comp_penalty = 0.0
        end

        phi_diff = eval_phi(candidate, candidate.point.mu) - eval_phi(it, it.point.mu)
        #phi_diff = eval_phi_diff(it, candidate, candidate.point.mu)
        return phi_diff + comp_penalty
    else
        return Inf
    end
end

function vector_product(LowerTriangularMat, vector::Vector)
    v1 = LowerTriangularMat * vector
    v2 = LowerTriangularMat' * vector
    return v1 + v2 - diag(LowerTriangularMat) .* vector
end

function hess_product(it::Class_iterate, vector::Vector)
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
    return s .* y + dir.y .* s * step_size + dir.s .* y * step_size - (it.point.mu + dir.mu * step_size)
end

function merit_function_predicted_reduction(it::Class_iterate, dir::Class_point, step_size::Float64)
    C_k = norm(comp(it),Inf)
    P_k = norm(comp_predicted(it, dir, step_size), Inf)
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
    return scaled_dual_feas(it, par) + norm(comp(it), Inf)
end
