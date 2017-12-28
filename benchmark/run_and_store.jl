function ipopt_solve(snlp::AbstractNLPModel, problem_name::String, test_name::String, my_par::OnePhase.Class_parameters)
    solver = IpoptSolver(print_level=3, tol=my_par.term.tol_opt, max_iter=3000, max_cpu_time=60.0^2, nlp_scaling_method="none", bound_relax_factor = 0.0, acceptable_iter=999999)

    println("RUNNING $problem_name")

    ORG_STDOUT = STDOUT
    file = open("../results/$(test_name)/log/$(problem_name).txt", "w")
    redirect_stdout(file)
    summary = problem_summary2()
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
      summary.status = status;
      summary.it_count = t;

      x = MathProgBase.getsolution(mp)
      set_cutest_info_ipopt!(summary, mp.inner, snlp, x)
    catch(e)
      println("Uncaught error in algorithm!!!")
      @show e;
      summary.status = :ERR
      summary.it_count = -1;
    end

    summary.total_time = time() - start_time;

    redirect_stdout(ORG_STDOUT)
    finalize(snlp)
    close(file)

    println("it count = ", summary.it_count)
    println("status = ", summary.status)

    save("../results/$(test_name)/summary.jld","summary",summary)

    return summary
end

function one_phase_solve(snlp::AbstractNLPModel, problem_name::String, test_name::String, my_par::OnePhase.Class_parameters)
    println("RUNNING $problem_name")
    ORG_STDOUT = STDOUT
    file = open("../results/$(test_name)/log/$(problem_name).txt", "w")
    redirect_stdout(file)
    summary = problem_summary2()
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

        @assert(is_feasible(intial_it, my_par.ls.comp_feas))
        iter, status, history, t, err = one_phase_IPM(intial_it, my_par, timer);

        pause_advanced_timer(timer)

        print_timer_stats(timer)

        #master_timer = merge_timers(timer, master_timer)

        save("../results/$(test_name)/jld/$(problem_name).jld","history",history, "timer", timer)

        set_info_me!(summary, history, status)
        #.it_count = t;
    catch(e)
        println("Uncaught error in algorithm!!!")
        @show e;

        if isa(e, Eval_NaN_error)
          summary.status = :NaN_ERR
          summary.it_count = -1;
        else
          summary.status = :ERR
          summary.it_count = -1;
        end
    end
    summary.total_time = time() - start_time;

    redirect_stdout(ORG_STDOUT)
    close(file)

    println("it count = ", summary.it_count)
    println("status = ", summary.status)

    timer_file = open("../results/$(test_name)/timer.txt", "w")
    #print_timer_stats(timer_file, master_timer)
    close(timer_file)

    return summary
end
