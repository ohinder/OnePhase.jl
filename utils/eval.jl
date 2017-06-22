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

function get_fval(it::Class_iterate)
    return it.cache.fval
end

function get_grad(it::Class_iterate)
    return it.cache.grad
end

function get_cons(it::Class_iterate)
    return it.cache.cons
end

function get_jac(it::Class_iterate)
    return it.cache.J
end

function get_primal_res(it::Class_iterate)
    return it.cache.cons - it.point.s
end

function get_max_vio(it::Class_iterate)
    return -min(0.0,minimum(it.cache.cons))
end

function eval_phi(it::Class_iterate)
    return it.cache.fval - it.point.mu * sum( log( it.point.s ) )
end

function eval_grad_phi(it::Class_iterate)
    y_tilde = it.point.mu ./ it.point.s
    return eval_grad_lag(it, y_tilde)
end


function get_lag_hess(it::Class_iterate)
    return it.cache.H
end


function eval_grad_lag(it::Class_iterate)
    return it.cache.grad - it.cache.J' * it.point.y
end

function eval_grad_lag(it::Class_iterate, y::Array{Float64,1})
    return it.cache.grad - it.cache.J' * y
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

        return eval_phi(it) + comp_penalty
    else
        return Inf
    end
end

function vector_product(LowerTriangularMat, vector)
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

    return step_size * dot(dir.x, eval_grad_phi(it)) + step_size^2 * 0.5 * (dot(dir.x, h_prod) + J_gain)
end

function comp_predicted(it::Class_iterate, dir::Class_point, step_size::Float64)
    y = it.point.y;
    s = it.point.s;
    return s .* y + dir.y .* s * step_size + dir.s .* y * step_size - it.point.mu
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
