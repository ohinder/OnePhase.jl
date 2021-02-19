println("GUROBI NOT WORKING")
#using Gurobi

function LP_init(nlp::Class_CUTEst, pars::Class_parameters, timer::class_advanced_timer)
    start_advanced_timer(timer, "INIT/x")
    x0 = suggested_starting_point(nlp)

    #@show norm(x0)

    if pars.init.start_satisfying_bounds
      x = projection_onto_bounds_ipopt_style( nlp, pars, x0 )
    end

    pause_advanced_timer(timer, "INIT/x")

    a, J, g = eval_init(nlp, pars, timer, x)

    s = a - 2 * minimum(a)
    @assert(all(s .>= 0.0))

    # make bounds exact.
    bis = bound_indicies(nlp)
    s[bis] = a[bis]
    @assert(all(s .>= 0.0))

    ncon = length(a)
    nvar = length(g)

    m = Model(solver=GurobiSolver())
    @variable(m, yVar[1:ncon] >= 0.0)
    @variable(m, muVar >= 0.0)

    y0 = 10.0 * ones(length(a)) # TEMP

    theta = 1e-4 # parameter
    penalty = theta
    scale = mean(abs.(J))
    sc_J = J / scale
    sc_g = g / scale
    #@show sc_J
    @show mean(abs.(J))

    if false
        @objective(m, Min, (sc_J' * yVar - sc_g)' * (sc_J' * yVar - sc_g) + penalty * (yVar - y0)' * (yVar - y0))
        status = solve(m)
    else
        @variable(m, z1[1:nvar] >= 0.0)
        @variable(m, z2[1:nvar] >= 0.0)
        @variable(m, t >= 0.0)
        @objective(m, Min, sum(z1[i]^2 for i = 1:nvar) + penalty * sum((yVar - y0).^2)) # + sum(log(yVar[j]) for j=1:ncon))

        if false
            for j = 1:ncon
                @constraint(m, yVar[j] <= muVar * 10.0 / s[j])
                @constraint(m, yVar[j] >= muVar * 0.1 / s[j])
            end
        end

        for i = 1:nvar
            @constraint(m, sum(sc_J[i,j] * yVar[j] for j = 1:nvar) - sc_g[i] <= z1[i])
            @constraint(m, sc_g[i] - sum(sc_J[i,j] * yVar[j] for j = 1:nvar) <= z1[i])
            #@constraint(m, sum(sc_J[i,j] * yVar[j] for j = 1:nvar) - sc_g[i] <= z2[i])
            #@constraint(m, sc_g[i] - sum(sc_J[i,j] * yVar[j] for j = 1:nvar) <= z2[i])
        end
        status = solve(m)
    end

    #mu = getvalue(muVar)
    y = getvalue(yVar)
    #@show getvalue(muVar)
    mu = - minimum(a) * mean(y) #/ length(y)
    #mu = mean(s .* y)
    obj = dot(J' * y - g, J' * y - g) + penalty * dot(y - y0, y - y0)
    @show maximum(y), getobjectivevalue(m), penalty, obj
    buffer = 2.0
    y_c = mu ./ s;
    @show norm(y_c,Inf), mu, minimum(s)
    y = min.( y_c / (pars.ls.comp_feas * buffer), max.(y, pars.ls.comp_feas * y_c * buffer)) # project onto complementarity constraints

    #@show mu, theta
    #error("")

    start_advanced_timer(timer, "INIT/construct_class")
    init_point = Class_point(x, y, s, mu)
    check_for_nan(init_point)
    init_it = Class_iterate(init_point, nlp, Class_local_info(), timer, pars);
    pause_advanced_timer(timer, "INIT/construct_class")

    return init_it
end
