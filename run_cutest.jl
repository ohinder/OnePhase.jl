include("include.jl")

#folder_name = ARGS[1]

function set_cutest_info_ipopt!(info::problem_summary, ipopt_solver, nlp_raw::CUTEst.CUTEstModel, x::Array{Float64})
  num_vars = length(nlp_raw.meta.lvar)
  x_true = x[1:num_vars]

  info.fval = obj(nlp_raw, x_true)
  a = cons(nlp_raw, x_true);
  info.con_vio = max(0.0, maximum(nlp_raw.meta.lvar - x_true), maximum(x_true - nlp_raw.meta.uvar), maximum(nlp_raw.meta.lcon - a), maximum(a - nlp_raw.meta.ucon))

  info.dual_feas = norm(grad(nlp_raw, x_true) + jac(nlp_raw, x_true)' * ipopt_solver.mult_g + ipopt_solver.mult_x_U - ipopt_solver.mult_x_L, Inf);
  info.comp = maximum(ipopt_solver.mult_x_U .* min(1e16, abs(nlp_raw.meta.uvar - x_true) ) ) + maximum( ipopt_solver.mult_x_L .* min(1e16, abs(x_true - nlp_raw.meta.lvar)) )
end

function run_cutest_problems_on_solver(problems::Array{String,1}, test_name::String, solver)
    summary = Dict{String, problem_summary}()

    if_mkdir("results/$test_name")
    if_mkdir("results/$test_name/log")

    for problem_name in problems
          println("RUNNING $problem_name")

          ORG_STDOUT = STDOUT
          file = open("results/$(test_name)/log/$(problem_name).txt", "w")
          redirect_stdout(file)
          nlp_raw = CUTEstModel(problem_name)
          summary[problem_name] = problem_summary()
          start_time = time()

          try
            t = 0
            function intermediate_ipopt(alg_mod::Int,iter_count::Int,obj_value::Float64,inf_pr::Float64, inf_du::Float64,  mu::Float64, d_norm::Float64, regularization_size::Float64, alpha_du::Float64, alpha_pr::Float64, ls_trials::Int)
              t = iter_count
              return true  # Keep going
            end

            mp = NLPModels.NLPtoMPB(nlp_raw, solver)
            setIntermediateCallback(mp.inner, intermediate_ipopt)
            MathProgBase.optimize!(mp)
            #x = MathProgBase.getsolution(mp)
            #solver = MathProgBase.getrawsolver(mp)
            status = MathProgBase.status(mp)
            summary[problem_name].status = status;
            summary[problem_name].it_count = t;

            x = MathProgBase.getsolution(mp)
            set_cutest_info_ipopt!(summary[problem_name], mp.inner, nlp_raw, x)
          catch(e)
            println("Uncaught error in algorithm!!!")
            @show e;
            summary[problem_name].status = :ERR
            summary[problem_name].it_count = -1;
          end

          summary[problem_name].total_time = time() - start_time;

          redirect_stdout(ORG_STDOUT)
          finalize(nlp_raw)
          close(file)

          println("it count = ", summary[problem_name].it_count)
          println("status = ", summary[problem_name].status)

          save("results/$(test_name)/summary.jld","summary",summary)

          summary_file = open("results/$(test_name)/summary.txt", "w")
          write_summary(summary_file, summary)
          close(summary_file)
    end
end


function run_cutest_problems_using_our_solver(problems::Array{String,1}, test_name::String, par::Class_parameters)

    if_mkdir("results/$test_name")
    if_mkdir("results/$test_name/log")
    if_mkdir("results/$test_name/jld")

    summary = Dict{String, problem_summary}()

    par_file = open("results/$(test_name)/par.txt", "w")
    write_pars(par_file, par)
    close(par_file)

    master_timer = class_advanced_timer()

    for problem_name in problems
          println("RUNNING $problem_name")
          ORG_STDOUT = STDOUT
          file = open("results/$(test_name)/log/$(problem_name).txt", "w")
          redirect_stdout(file)
          nlp_raw = CUTEstModel(problem_name)
          summary[problem_name] = problem_summary()
          start_time = time()

          try
              nlp = Class_CUTEst(nlp_raw)

              timer = class_advanced_timer()
              start_advanced_timer(timer)
              #include("include.jl")
              #intial_it = initial_point_satisfy_bounds(nlp, my_par)
              start_advanced_timer(timer, "INIT")
              intial_it = init(nlp, my_par, timer)
              pause_advanced_timer(timer, "INIT")

              #intial_it = initial_point_generic(nlp, my_par, nlp_raw.meta.x0)

              @assert(is_feasible(intial_it, my_par.comp_feas))
              iter, status, history, t, err = one_phase_IPM(intial_it, my_par, timer);

              pause_advanced_timer(timer)

              print_timer_stats(timer)

              master_timer = merge_timers(timer, master_timer)

              save("results/$(test_name)/jld/$(problem_name).jld","history",history, "timer", timer)

              summary[problem_name].status = status;
              set_info_me!(summary[problem_name], history)
              #.it_count = t;
          catch(e)
              println("Uncaught error in algorithm!!!")
              @show e;
              summary[problem_name].status = :ERR
              summary[problem_name].it_count = -1;
          end
          summary[problem_name].total_time = time() - start_time;

          redirect_stdout(ORG_STDOUT)
          finalize(nlp_raw)
          close(file)

          println("it count = ", summary[problem_name].it_count)
          println("status = ", summary[problem_name].status)

          save("results/$(test_name)/summary.jld","summary",summary, "pars", par)

          summary_file = open("results/$(test_name)/summary.txt", "w")
          write_summary(summary_file, summary)
          close(summary_file)

          timer_file = open("results/$(test_name)/timer.txt", "w")
          print_timer_stats(timer_file, master_timer)
          close(timer_file)
    end
end

#function run_cutest_problems_using_IPOPT(problems::Array{String,1}, test_name::String)
#
#end

#run_cutest_problems(["DISCS"], "test")


function filter_cutest(problem)
    # medium
    correct_size = 50 <= problem["variables"]["number"] + problem["constraints"]["number"] && problem["constraints"]["number"] >= 10 && problem["variables"]["number"] <= 600 && problem["constraints"]["number"] <= 1000
    regular = problem["derivative_order"] >= 2 && problem["regular"] == true
    # small
    #correct_size = 100 <= problem["variables"]["number"] + problem["constraints"]["number"] && problem["constraints"]["number"] >= 10 && problem["variables"]["number"] + problem["constraints"]["number"] <= 300
    #regular = problem["derivative_order"] >= 2 && problem["regular"] == true
    if correct_size && regular
        return true
    else
      return false
    end
end

#problem_list = CUTEst.select(max_var=100, max_var=1000, min_con=100, max_con=3000)
problem_list = CUTEst.select(custom_filter=filter_cutest)
problem_list = convert(Array{String,1},problem_list)

if false
problem_list = ["PT", "AGG"]
folder_name = "test_run"
if_mkdir("results/$folder_name")
run_cutest_problems_using_our_solver(problem_list, "$folder_name", my_par)
end

if true
    folder_name = "mehotra_intial_point6"
    if_mkdir("results/$folder_name")
    run_cutest_problems_using_our_solver(problem_list, "$folder_name", my_par)

    #folder_name = "par2/mehotra-no-satisfy"
    #if_mkdir("results/$folder_name")
    #my_par.start_satisfying_bounds = false
    #run_cutest_problems_using_our_solver(problem_list, "$folder_name", my_par)
end

if false
    folder_name = "par_hess1"
    if_mkdir("results/$folder_name")

    for mu_ratio in [0.01, 1.0, 100.0]
        my_par.mu_primal_ratio = mu_ratio
        run_cutest_problems_using_our_solver(problem_list, "$folder_name/mu_ratio-$mu_ratio", my_par)
    end
end

if false
  using Ipopt
  my_solver = IpoptSolver(print_level=3)
  run_cutest_problems_on_solver(problem_list, "ipopt_test2", my_solver)
end
