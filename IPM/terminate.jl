function unboundedness_measure(iter::Class_iterate, tol::Float64)
    return norm(iter.point.x) #max( tol,get_max_vio(iter) ) / (tol * max(1.0,-get_fval(iter)))
end

function terminate(iter::Class_iterate, par::Class_parameters)
    y = iter.point.y
    fark_feas1 = farkas_certificate(iter, y)
    fark_feas2 = stationary_infeasible_measure(iter, y)
    unboun = unboundedness_measure(iter, par.term.tol_opt)

    max_vio = get_max_vio(iter)

    comp_scaled = maximum(iter.point.s .* y) * dual_scale(iter, par)

    if scaled_dual_feas(iter, par) < par.term.tol_opt && comp_scaled < par.term.tol_opt && max_vio < par.term.tol_opt
        return :optimal
    elseif max_vio > par.term.tol_opt && fark_feas1 < par.term.tol_inf_1 && fark_feas2 < par.term.tol_inf_2  #&& norm(get_y(iter), Inf) > 1/tol #&& norm(get_y(iter),Inf) > 1/tol
        return :primal_infeasible
    elseif unboun < par.term.tol_unbounded # TMP!!!
        return :dual_infeasible
    elseif norm(get_grad(iter),1) > par.term.max_gradient
        return :max_gradient
    else
        return false
    end
end
