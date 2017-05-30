include("include.jl")

folder_name = ARGS[1]

function run_cutest_problems_using_our_solver(problems::Array{String,1}, test_name::String, par::Class_parameters)

    if_mkdir("results/$test_name")
    if_mkdir("results/$test_name/log")
    if_mkdir("results/$test_name/jld")

    summary = Dict{String, problem_summary}()

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

              reset_advanced_timer()
              start_advanced_timer()
              intial_it = init(nlp, par)
              iter, status, history, t = one_phase_IPM(intial_it, par);
              pause_advanced_timer()
              print_timer_stats()

              save("results/$(test_name)/jld/$(problem_name).jld","history",history)

              summary[problem_name].status = status;
              summary[problem_name].it_count = t;
          catch(e)
              @show e;
              summary[problem_name].status = :ERROR
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

function run_cutest_problems_using_IPOPT(problems::Array{String,1}, test_name::String)

end

#run_cutest_problems(["DISCS"], "test")

function filter_cutest(problem)
    correct_size = 50 <= problem["variables"]["number"] + problem["constraints"]["number"] && problem["constraints"]["number"] >= 10 && problem["variables"]["number"] <= 600 && problem["constraints"]["number"] <= 1000
    regular = problem["derivative_order"] >= 2 && problem["regular"] == true
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
problem_list = ["PT"]
folder_name = "test_run"
end

if true
    run_cutest_problems_using_our_solver(problem_list, "$folder_name", my_par)
end

if false
    folder_name = "par1"
    if_mkdir("results/$folder_name")

    for mu_ratio in [0.01, 1.0, 100.0]
        my_par.mu_primal_ratio = mu_ratio
        run_cutest_problems_using_our_solver(problem_list, "$folder_name/mu_ratio-$mu_ratio", my_par)
    end

    my_par.mu_primal_ratio = 1.0
    my_par.start_satisfying_bounds = true
    run_cutest_problems_using_our_solver(problem_list, "$folder_name/bounds", my_par)

    my_par.start_satisfying_bounds = false
    run_cutest_problems_using_our_solver(problem_list, "$folder_name/nobounds", my_par)
end
