function comp(it::Class_iterate)
    return it.point.s .* it.point.y - it.point.mu
end

function eval_f(it::Class_iterate)
    start_advanced_timer("eval/f")
    fval = eval_f(it.nlp, it.point.x)
    pause_advanced_timer("eval/f")
    return fval
end

function eval_a(it::Class_iterate)
    start_advanced_timer("eval/a")
    a = eval_a(it.nlp, it.point.x)
    pause_advanced_timer("eval/a")
    return a
end

function eval_primal_residual(it)
    primal_r = eval_a(it) - it.point.s
    return primal_r
end

function eval_grad_f(it::Class_iterate)
    start_advanced_timer("eval/grad")
    grad = eval_grad_lag(it.nlp, it.point.x, zeros(ncon(it)), 1.0)
    pause_advanced_timer("eval/grad")
    return grad
end

function eval_phi(it::Class_iterate)
    phi = eval_f(it) - it.point.mu * sum(log(it.point.s))
    return phi
end

function eval_grad_phi(it::Class_iterate)
    start_advanced_timer("eval/grad_phi")
    grad_phi = eval_grad_lag(it.nlp, it.point.x, it.point.mu ./ it.point.s, 1.0)
    pause_advanced_timer("eval/grad_phi")
    return grad_phi
end

function eval_grad_lag(it::Class_iterate, w::Float64=1.0)
    start_advanced_timer("eval/grad_lag")
    grad_lag = eval_grad_lag(it.nlp, it.point.x, it.point.y, w)
    pause_advanced_timer("eval/grad_lag")
    return grad_lag
end

function eval_farkas(it::Class_iterate)
    eval_farkas = eval_grad_lag(it,0.0)/norm(it.point.y,Inf)
    return eval_farkas
end

function eval_jac(it::Class_iterate)
    start_advanced_timer("eval/jac")
    jac = eval_jac(it.nlp, it.point.x)
    pause_advanced_timer("eval/jac")

    return jac
end


function eval_merit_function(it::Class_iterate)
    merit = eval_phi(it) + norm(comp(it),Inf)^3 / (it.point.mu)^2

    return merit
end

function eval_lag_hess(it::Class_iterate, w::Float64=1.0)
    start_advanced_timer("eval/hess")
    lag_hess = eval_lag_hess(it.nlp, it.point.x, it.point.y, w)
    pause_advanced_timer("eval/hess")

    return lag_hess
end

function comp_predicted(it::Class_iterate, dir::Class_point)
    y = it.point.y;
    s = it.point.s;
    return s .* y + dir.y .* s + dir.s .* y - it.point.mu
end

function phi_predicted_reduction(it::Class_iterate, dir::Class_point)
    grad_phi = eval_grad_phi(it)
    return dot(dir.x,grad_phi)
end

function merit_function_predicted_reduction(it::Class_iterate, dir::Class_point)
    C_k = norm(comp(it),Inf)
    P_k = norm(comp_predicted(it, dir), Inf)
    #@show C_k, P_k
    comp_penalty = (P_k^3 - C_k^3) / (it.point.mu)^2
    #@show comp_penalty
    predict_red = phi_predicted_reduction(it, dir)
    #@show predict_red

    return predict_red + comp_penalty
end
