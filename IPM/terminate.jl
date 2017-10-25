function terminate(iter::Class_iterate, par::Class_parameters)
    tol = par.tol
    #a_neg_part = max(-eval_a(iter),0.0)
    #y = get_y(iter)
    #J = eval_jac(iter)
    fark_feas1 = eval_farkas(iter)
    #fark_feas1 = norm(y' * J, Inf) / min(norm(y, Inf), maximum(a_neg_part .* y))
    #fark_feas2 = norm(eval_grad_lag(iter, 0.0), Inf) /  min(norm(y, Inf), maximum(a_neg_part .* y))
    #fark_feas3 = norm(eval_grad_lag(iter, 1.0), Inf) /  min(norm(y, Inf), maximum(a_neg_part .* y))

    #@show fark_feas1, fark_feas2, fark_feas3
    max_vio = get_max_vio(iter)

    comp_scaled = maximum(iter.point.s .* iter.point.y) * dual_scale(iter, par)

    if scaled_dual_feas(iter, par) < tol && comp_scaled < tol && max_vio < tol
        return :optimal
    elseif fark_feas1 < tol && max_vio > tol #&& norm(get_y(iter), Inf) > 1/tol #&& norm(get_y(iter),Inf) > 1/tol
        return :primal_infeasible
    elseif max( 1.0, max_vio )/ min(max(1.0,-get_fval(iter)),norm(get_x(iter),Inf)) < tol # TMP!!!
        return :dual_infeasible
    else
        return false
    end
end
