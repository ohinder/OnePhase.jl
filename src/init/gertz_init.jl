function gertz_init(nlp::Class_CUTEst, pars::Class_parameters, timer::class_advanced_timer)
    #println("0'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    start_advanced_timer(timer, "INIT/x")
    #println("1'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    x0 = suggested_starting_point(nlp)
    #println("2'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    #x0 += rand(length(x0))

    if pars.init.start_satisfying_bounds
      #println("3'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
      x = projection_onto_bounds_ipopt_style( nlp, pars, x0 )
    end
    pause_advanced_timer(timer, "INIT/x")
    #println("4'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    y = ones(ncon(nlp))
    #println("5'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    a, J, g = eval_init(nlp, pars, timer, x)
    #println("6'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    #s = max(a, ones(ncon(nlp)))
    s_thres = 1e-4
    #println("1''HHHHHHHHHHHHH,$s_thres,HHHHHHHHHHH,$a,HHHHHHHHHH,$minimum(a)")
    #s_thres = LinearAlgebra.norm(g - J' * y, Inf)
    #s = max(a,s_thres) #
    d_s = max(s_thres,-2.0 * minimum(a))
    #println("2''HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    s = a .+ d_s
    #println("3''HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    mu = max(s_thres,-2.0 * minimum(a))
    #y = ones(ncon(nlp)) / d_s
    #y = mu ./ s
    #println("4''HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    init_point = Class_point(x, y, s, mu)
    #println("7'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    iter = Class_iterate(init_point, nlp, Class_local_info(), timer, pars);
    #println("8'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    kkt_solver = pick_KKT_solver(pars);
    #println("9'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    initialize!(kkt_solver, iter)
    #println("10'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    form_system!(kkt_solver, iter, timer)
    #println("11'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    fact_succeed, inertia_num_fac, new_delta = ipopt_strategy!(iter, kkt_solver, pars, timer)
    #println("12'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    kkt_associate_rhs!(kkt_solver, iter, Reduct_affine(), timer)
    #println("13'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    compute_direction!(kkt_solver, timer)
    #println("14'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    dir = kkt_solver.dir
    y_temp = init_point.y + dir.y
    s_temp = -a
    #println("0HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    s, y = mehortra_guarding( nlp, pars, timer, x, y_temp, s_temp, a, J, g )
    #println("1HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    @assert(all(s .>= 0.0))
    #println("2HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    @assert(all(s .> 0.0))
    #println("3HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    @assert(all(y .> 0.0))
    #println("4HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    mu = Statistics.mean(s .* y)
    @assert(mu > 0.0)
    #println("5HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    conWeight = ((s - a) / mu)
    @assert(all(conWeight .>= 0.0))
    #println("6HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    @assert(all(a + mu * conWeight .>= 0.0))
    #println("7HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
    return correct_guess3( nlp, pars, timer, x, a, J, g, y, mu, conWeight )
end

function center_dual!(iter::Class_iterate, pars::Class_parameters)
    y_c = iter.point.mu ./ iter.point.s
    y = iter.point.y
    y = min.( y_c / pars.ls.comp_feas_agg, max.(y, pars.ls.comp_feas_agg * y_c)) # project onto complementarity constraints
    iter.point.y = y;
end
