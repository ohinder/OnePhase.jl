function gertz_init(nlp::Class_CUTEst, pars::Class_parameters, timer::class_advanced_timer)
    start_advanced_timer(timer, "INIT/x")
    x0 = suggested_starting_point(nlp)
    if pars.init.start_satisfying_bounds
      x = projection_onto_bounds_ipopt_style( nlp, pars, x0 )
    end
    pause_advanced_timer(timer, "INIT/x")
    y = ones(ncon(nlp))
    a, J, g = eval_init(nlp, pars, timer, x)

    #s = max(a, ones(ncon(nlp)))
    s_thres = 1e-4
    #s_thres = LinearAlgebra.norm(g - J' * y, Inf)
    #s = max(a,s_thres) #
    d_s = max(s_thres,-2.0 * minimum(a))
    s = a .+ d_s
    mu = max(s_thres,-2.0 * minimum(a))
    #y = ones(ncon(nlp)) / d_s
    #y = mu ./ s
    init_point = Class_point(x, y, s, mu)
    iter = Class_iterate(init_point, nlp, Class_local_info(), timer, pars);
    kkt_solver = pick_KKT_solver(pars);
    initialize!(kkt_solver, iter)
    form_system!(kkt_solver, iter, timer)
    fact_succeed, inertia_num_fac, new_delta = ipopt_strategy!(iter, kkt_solver, pars, timer)
    kkt_associate_rhs!(kkt_solver, iter, Reduct_affine(), timer)
    compute_direction!(kkt_solver, timer)
    dir = kkt_solver.dir
    y_temp = init_point.y + dir.y
    s_temp = -a
    s, y = mehortra_guarding(nlp, pars, timer, x, y_temp, s_temp, a, J, g)
    @assert(all(s .>= 0.0))
    @assert(all(s .> 0.0))
    @assert(all(y .> 0.0))
    mu = Statistics.mean(s .* y)
    @assert(mu > 0.0)
    conWeight = ((s - a) / mu)
    @assert(all(conWeight .>= 0.0))
    @assert(all(a + mu * conWeight .>= 0.0))
    return correct_guess3( nlp, pars, timer, x, a, J, g, y, mu, conWeight )
end

function center_dual!(iter::Class_iterate, pars::Class_parameters)
    y_c = iter.point.mu ./ iter.point.s
    y = iter.point.y
    y = min.( y_c / pars.ls.comp_feas_agg, max.(y, pars.ls.comp_feas_agg * y_c)) # project onto complementarity constraints
    iter.point.y = y;
end
