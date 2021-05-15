

function terminate(iter::Class_iterate, par::Class_parameters)
    y = iter.point.y
    fark_feas1 = farkas_certificate(iter, y)
    fark_feas2 = stationary_infeasible_measure(iter, y)

    max_vio = get_max_vio(iter) # ???

    comp_scaled = maximum(iter.point.s .* y) * dual_scale(iter, par)
    
    if scaled_dual_feas(iter, 0.0, par) < par.term.tol_opt && comp_scaled < par.term.tol_opt && max_vio < par.term.tol_opt
        return :Optimal
    elseif max_vio > par.term.tol_opt && fark_feas1 < par.term.tol_inf_1 && fark_feas2 < par.term.tol_inf_2  #&& LinearAlgebra.norm(get_y(iter), Inf) > 1/tol #&& LinearAlgebra.norm(get_y(iter),Inf) > 1/tol
        return :primal_infeasible
    elseif LinearAlgebra.norm(iter.point.x,Inf) > 1.0 / par.term.tol_unbounded # TMP!!!
        return :dual_infeasible
    elseif LinearAlgebra.norm(iter.cache.grad,Inf) > par.term.grad_max
        return :max_gradient
    else
        return false
    end
end
