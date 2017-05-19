function test(it::Class_iterate)
    @show eval_f(it)
    @show eval_phi(it)
    @show eval_a(it)
    @show comp(it)
    @show eval_merit_function(it)

    @show eval_grad_lag(it)
    @show eval_grad_phi(it)

    @show eval_lag_hess(it)
end

function run_lp()
    #intial_it = toy_LP()
    intial_it = toy_LP2()

    test(intial_it);
    validate(intial_it);

    dir, = compute_direction(intial_it, intial_it, speye(length(intial_it.point.x)), par, :affine)
    shrink_direction!(dir,0.1)

    is_feasible(intial_it, par.comp_feas)

    dir = zero_point(length(intial_it.point.x),length(intial_it.point.s))
    acceptable(intial_it, intial_it, dir, par)

    @show eval_primal_residual(intial_it)
    final_iter, status = one_phase(intial_it, par)
    @show eval_primal_residual(final_iter)

    simple_LP_solver(intial_it, par)
end

function initial_point(nlp::abstract_nlp, x::Array{Float64,1})
    a = eval_a(nlp, x)
    infeas = norm(min(a,0.0),Inf)
    s = max(0, a) + infeas
    mu = maximum(s - a)
    y = mu ./ s;



    init_point = Class_point(x, y, s, mu)

    initial_it = Class_iterate(init_point, nlp, Class_local_info(0.0));
    #A = eval_jac(initial_it)
    #grad = eval_
    return initial_it
end

function comp_matrix(n::Int64)
    M = spzeros(n,n)
    for i = 1:(n-1)
        M[i,i+1] = 1.0
        M[i+1,i] = 1.0
    end
    return M
end

function run_netlib_lp(name::String, par::Class_parameters)
    LP = read_lp(name)

    #LP.Q = 1e-1 * comp_matrix(dim(LP))
    #LP.b += 1e-5

    println("Running:", name)
    @show dim(LP), ncon(LP)

    x = zeros(dim(LP));

    intial_it = initial_point(LP, x)

    iter, status, hist = one_phase_IPM(intial_it, par)
end
