include("run_cutest.jl")

function perturb_cons(nr::CUTEstModel,scale::Float64)
    srand(1)

    lcon = deepcopy(nr.meta.lcon);
    ucon = deepcopy(nr.meta.ucon);
    lvar = deepcopy(nr.meta.lvar);
    uvar = deepcopy(nr.meta.uvar);

    f(x) = obj(nr, x)
    #shift = scale * (rand(length(lcon)) - 0.5)
    shift = scale * ones(length(lcon)) #(rand(length(lcon)) - 0.5)
    c(x) = cons(nr, x) #
    g(x) = grad(nr, x)
    J(x) = jac(nr, x)
    H(x; obj_weight=1.0, y=zeros) = hess(nr, x, obj_weight=obj_weight, y=y)
    Hcoord(x; obj_weight=1.0, y=zeros) = hess_coord(nr, x, obj_weight=obj_weight, y=y)
    Jcoord(x) = jac_coord(nr, x)

    snlp = SimpleNLPModel(f, nr.meta.x0;
    lvar = lvar, uvar = uvar,
    lcon = lcon - shift, ucon = ucon - shift,
    name = "Generic",
    g=g,
    c=c,
    J=J,
    H=H,
    Jcoord=Jcoord,
    Hcoord=Hcoord)

    return snlp
end

function run_infeas(test_name, problem_list, my_par, solve_func!)
    summary = Dict{String, problem_summary2}()

    if_mkdir("../results/$test_name")
    if_mkdir("../results/$test_name/log")
    if_mkdir("../results/$test_name/jld")

    for problem_name in problem_list
        nr = CUTEstModel(problem_name)
        snlp = perturb_cons(nr, 1.0)
        summary[problem_name] = solve_func!(snlp, problem_name, test_name, my_par)
        finalize(nr)

        write_summary("../results/$(test_name)/summary.txt", summary)
        save("../results/$(test_name)/summary.jld","summary",summary, "pars", my_par)
    end
end

problem_list = get_problem_list(50,10000)

my_par = Class_parameters()
run_infeas("one_phase/infeas-Dec20", problem_list, my_par, one_phase_solve)
#run_infeas("ipopt/infeas-test", problem_list, my_par, ipopt_solve)
