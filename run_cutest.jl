include("include.jl")

test_name = ARGS[1]

type problem_summary
    status::Symbol
    it_count::Int64
    total_time::Float64

    function problem_summary()
        return new()
    end
end

function write_summary(stream::IOStream, summary::Dict{String, problem_summary})
    write(stream, "problem_name \t status \t it count \t total_time \n");
    for (problem_name, info) in summary
      write(stream, "$(pd(problem_name)) \t $(pd(info.status)) \t $(pd(info.it_count)) \t $(rd(info.total_time))\n");
    end
end



function run_cutest_problems(problems::Array{String,1}, test_name::String)

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
              intial_it = initial_point_satisfy_bounds(nlp, nlp_raw.meta.x0, my_par)

              iter, status, history, t = one_phase_IPM(intial_it, my_par);
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

#run_cutest_problems(["DISCS"], "test")

function filter_cutest(problem)
    if problem["derivative_order"] >= 2 && 100 <= problem["variables"]["number"] && problem["variables"]["number"] <= 1000 && 100 <= problem["constraints"]["number"] && problem["constraints"]["number"] <= 3000
        return true
    else
      return false
    end
end

#problem_list = CUTEst.select(max_var=100, max_var=1000, min_con=100, max_con=3000)
problem_list = CUTEst.select(custom_filter=filter_cutest)
problem_list = convert(Array{String,1},problem_list)
#problem_list = ["PT"] #"DISCS"]
run_cutest_problems(problem_list, test_name)
