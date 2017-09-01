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
    c(x) = cons(nr, x) - shift
    g(x) = grad(nr, x)
    J(x) = jac(nr, x)
    H(x; obj_weight=1.0, y=zeros) = hess(nr, x, obj_weight=obj_weight, y=y)
    Hcoord(x; obj_weight=1.0, y=zeros) = hess_coord(nr, x, obj_weight=obj_weight, y=y)
    Jcoord(x) = jac_coord(nr, x)

    snlp = SimpleNLPModel(f, nr.meta.x0;
    lvar = lvar, uvar = uvar,
    lcon = lcon, ucon = ucon,
    name = "Generic",
    g=g,
    c=c,
    J=J,
    H=H,
    Jcoord=Jcoord,
    Hcoord=Hcoord)

    return snlp
end

function ipopt_solve!(snlp::AbstractNLPModel, problem_name, test_name, my_par, summary)
    solver = IpoptSolver(print_level=3, tol=1e-6, max_iter=3000, max_cpu_time=60.0, nlp_scaling_method="none", bound_relax_factor = 0.0, acceptable_iter=999999)

    println("RUNNING $problem_name")

    ORG_STDOUT = STDOUT
    file = open("../results/$(test_name)/log/$(problem_name).txt", "w")
    redirect_stdout(file)
    summary[problem_name] = problem_summary()
    start_time = time()

    try
      t = 0
      function intermediate_ipopt(alg_mod::Int,iter_count::Int,obj_value::Float64,inf_pr::Float64, inf_du::Float64,  mu::Float64, d_norm::Float64, regularization_size::Float64, alpha_du::Float64, alpha_pr::Float64, ls_trials::Int)
        t = iter_count
        return true  # Keep going
      end

      mp = NLPModels.NLPtoMPB(snlp, solver)
      setIntermediateCallback(mp.inner, intermediate_ipopt)
      MathProgBase.optimize!(mp)
      #x = MathProgBase.getsolution(mp)
      #solver = MathProgBase.getrawsolver(mp)
      status = MathProgBase.status(mp)
      summary[problem_name].status = status;
      summary[problem_name].it_count = t;

      x = MathProgBase.getsolution(mp)
      set_cutest_info_ipopt!(summary[problem_name], mp.inner, snlp, x)
    catch(e)
      println("Uncaught error in algorithm!!!")
      @show e;
      summary[problem_name].status = :ERR
      summary[problem_name].it_count = -1;
    end

    summary[problem_name].total_time = time() - start_time;

    redirect_stdout(ORG_STDOUT)
    finalize(snlp)
    close(file)

    println("it count = ", summary[problem_name].it_count)
    println("status = ", summary[problem_name].status)

    save("../results/$(test_name)/summary.jld","summary",summary)

    summary_file = open("../results/$(test_name)/summary.txt", "w")
    write_summary(summary_file, summary)
    close(summary_file)
end


#=
if_mkdir("../results/$test_name")
if_mkdir("../results/$test_name/log")
if_mkdir("../results/$test_name/jld")

if isfile("../results/$(test_name)/summary.jld")
  summary = load("../results/$(test_name)/summary.jld","summary")
else
  summary = Dict{String, problem_summary}()
end

already_solved_problems = keys(summary)

par_file = open("../results/$(test_name)/par.txt", "w")
write_pars(par_file, par)
close(par_file)=#

function one_phase_solve!(snlp::AbstractNLPModel, problem_name, test_name, my_par, summary)
    println("RUNNING $problem_name")
    ORG_STDOUT = STDOUT
    file = open("../results/$(test_name)/log/$(problem_name).txt", "w")
    redirect_stdout(file)
    summary[problem_name] = problem_summary()
    start_time = time()

    try
        nlp = Class_CUTEst(snlp)

        timer = class_advanced_timer()
        start_advanced_timer(timer)
        #include("include.jl")
        #intial_it = initial_point_satisfy_bounds(nlp, my_par)
        start_advanced_timer(timer, "INIT")
        intial_it = init(nlp, my_par, timer)
        pause_advanced_timer(timer, "INIT")

        @assert(is_feasible(intial_it, my_par.comp_feas))
        iter, status, history, t, err = one_phase_IPM(intial_it, my_par, timer);

        pause_advanced_timer(timer)

        print_timer_stats(timer)

        #master_timer = merge_timers(timer, master_timer)

        save("../results/$(test_name)/jld/$(problem_name).jld","history",history, "timer", timer)

        summary[problem_name].status = status;
        set_info_me!(summary[problem_name], history)
        #.it_count = t;
    catch(e)
        println("Uncaught error in algorithm!!!")
        @show e;

        if isa(e, Eval_NaN_error)
          summary[problem_name].status = :NaN_ERR
          summary[problem_name].it_count = -1;
        else
          summary[problem_name].status = :ERR
          summary[problem_name].it_count = -1;
        end
    end
    summary[problem_name].total_time = time() - start_time;

    redirect_stdout(ORG_STDOUT)
    close(file)

    println("it count = ", summary[problem_name].it_count)
    println("status = ", summary[problem_name].status)

    save("../results/$(test_name)/summary.jld","summary",summary, "pars", my_par)

    summary_file = open("../results/$(test_name)/summary.txt", "w")
    write_summary(summary_file, summary)
    close(summary_file)

    timer_file = open("../results/$(test_name)/timer.txt", "w")
    #print_timer_stats(timer_file, master_timer)
    close(timer_file)
end



using Ipopt

function run_infeas(test_name, problem_list, my_par, solve_func!)
    summary = Dict{String, problem_summary}()

    if_mkdir("../results/$test_name")
    if_mkdir("../results/$test_name/log")
    if_mkdir("../results/$test_name/jld")

    for problem_name in problem_list
        nr = CUTEstModel(problem_name)
        snlp = perturb_cons(nr, 1.0)
        solve_func!(snlp, problem_name, test_name, my_par, summary)
        finalize(nr)
    end
end

problem_list = get_problem_list(100,1000)

run_infeas("ipopt/infeas-3", problem_list, my_par, ipopt_solve!)
run_infeas("one_phase/infeas-3", problem_list, my_par, one_phase_solve!)
