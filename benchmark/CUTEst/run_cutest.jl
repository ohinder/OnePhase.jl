include("../benchmark.jl")
using OnePhase, advanced_timer

function run_cutest_problems_on_solver(problems::Array{String,1}, test_name::String, solver)
    summary = Dict{String, problem_summary2}()

    if_mkdir("../results/$test_name")
    if_mkdir("../results/$test_name/log")

    for problem_name in problems
          println("RUNNING $problem_name")

          ORG_STDOUT = STDOUT
          file = open("../results/$(test_name)/log/$(problem_name).txt", "w")
          redirect_stdout(file)
          nlp_raw = CUTEstModel(problem_name)
          summary[problem_name] = problem_summary2()
          tic()

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
            println("here4")
          catch(e)
            println("Uncaught error in algorithm!!!")
            @show e;
            summary[problem_name].status = :ERR
            summary[problem_name].it_count = -1;
          end

          summary[problem_name].total_time = toc();

          redirect_stdout(ORG_STDOUT)
          finalize(nlp_raw)
          close(file)

          println("it count = ", summary[problem_name].it_count)
          println("status = ", summary[problem_name].status)

          save("../results/$(test_name)/summary.jld","summary",summary)

          summary_file = open("../results/$(test_name)/summary.txt", "w")
          write_summary(summary_file, summary)
          close(summary_file)
    end
end


function run_cutest_problems_using_our_solver(problems::Array{String,1}, test_name::String, par::Class_parameters)

    if_mkdir("../results/$test_name")
    if_mkdir("../results/$test_name/log")
    if_mkdir("../results/$test_name/jld")

    if isfile("../results/$(test_name)/summary.jld")
      summary = load("../results/$(test_name)/summary.jld","summary")
    else
      summary = Dict{String, problem_summary2}()
    end

    already_solved_problems = keys(summary)

    par_file = open("../results/$(test_name)/par.txt", "w")
    write_pars(par_file, par)
    close(par_file)

    master_timer = class_advanced_timer()

    for problem_name in problems
      if problem_name in already_solved_problems
          println("$problem_name already solved")
      else
          println("RUNNING $problem_name")
          ORG_STDOUT = STDOUT
          file = open("../results/$(test_name)/log/$(problem_name).txt", "w")
          redirect_stdout(file)
          summary[problem_name] = problem_summary2()

          nlp_raw = false
          nlp = false

          try
              try
                nlp_raw = CUTEstModel(problem_name)
                #nlp = Class_CUTEst(nlp_raw)
              catch(e)
                println("$problem_name failed to load")
                throw(e)
              end

              tic()

              #=timer = class_advanced_timer()
              start_advanced_timer(timer)
              #include("include.jl")
              #intial_it = initial_point_satisfy_bounds(nlp, par)
              start_advanced_timer(timer, "INIT")
              intial_it = OnePhase.init(nlp, par, timer)
              pause_advanced_timer(timer, "INIT")

              #intial_it = initial_point_generic(nlp, par, nlp_raw.meta.x0)

              @assert(OnePhase.is_feasible(intial_it, par.ls.comp_feas))
              iter, status, history, t, err = OnePhase.one_phase_IPM(intial_it, par, timer);

              pause_advanced_timer(timer)

              print_timer_stats(timer)
              =#
              iter, status, history, t, err, timer = OnePhase.one_phase_solve(nlp_raw, par)

              master_timer = merge_timers(timer, master_timer)

              save("../results/$(test_name)/jld/$(problem_name).jld","history",history, "timer", timer)

              set_info_me!(summary[problem_name], history, status)
              #.it_count = t;

              summary[problem_name].total_time = toc();
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

          redirect_stdout(ORG_STDOUT)
          finalize(nlp_raw)
          close(file)

          println("it count = ", summary[problem_name].it_count)
          println("status = ", summary[problem_name].status)

          save("../results/$(test_name)/summary.jld","summary",summary, "pars", par)

          summary_file = open("../results/$(test_name)/summary.txt", "w")
          write_summary(summary_file, summary)
          close(summary_file)

          timer_file = open("../results/$(test_name)/timer.txt", "w")
          print_timer_stats(timer_file, master_timer)
          close(timer_file)
      end
    end
end

#function run_cutest_problems_using_IPOPT(problems::Array{String,1}, test_name::String)
#
#end

#run_cutest_problems(["DISCS"], "test")






function test_problems(problem_list::Array{String,1}, start::Int64)
    for i = start:length(problem_list)
      problem_name = problem_list[i]
      println(problem_name)
      nlp_raw = CUTEstModel(problem_name)
      println(i)
      finalize(nlp_raw)
    end
end

function filter_string_array(problem_list::Array{String,1}, remove_list::Array{String,1})
    new_problem_list = Array{String,1}()
    for problem_name in problem_list
      if !(problem_name in remove_list)
        push!(new_problem_list, problem_name)
      end
    end
    return new_problem_list
end


function get_problem_list(filter::Function)
  #problem_list = CUTEst.select(max_var=100, max_var=1000, min_con=100, max_con=3000)
  problem_list = CUTEst.select(custom_filter=filter);
  problem_list = convert(Array{String,1},problem_list);

  remove_list = ["MISRA1D","OSBORNE1","LANCZOS2","MEYER3NE","ROSZMAN1","INTEQNE","MGH10",
  "ECKERLE4","RAT43","MISRA1A","KOWOSBNE","MOREBVNE",
  "JENSMPNE","MISRA1C","VARDIMNE","SANTA","BARDNE","MGH10S",
  "THURBER","MGH17S","BOX3NE","PENLT1NE","GAUSS2","GULFNE","BA-L1SP","ENSO","CHWIRUT1", "CHWIRUT2","NELSON","HAHN1","KIRBY2","MUONSINE","OSBORNE2","BENNETT5","BA-L1","GAUSS3",
  "GAUSS1","LANCZOS3","BIGGS6NE","LANCZOS1","VANDANIUMS","MISRA1B","MGH09","MGH17","WATSONNE",
  "DMN37142LS","VESUVIOU","PENLT2NE", "DMN37143LS","INTEQNELS","BROYDN3DLS","DMN15332LS",
  "DMN15102LS","BROYDNBDLS","VESUVIA","BA-L1SPLS","BA-L1SPLS","SANTALS","VESUVIO","DMN15333LS",
  "ARGTRIGLS","RAT42","SSINE","LSC1","LSC2","BOXBOD","POWELLSE","FREURONE","DANWOOD","HELIXNE",
  "DMN37142","DMN37143","DMN15102","DMN15333","DMN15332","DMN15103",
  "VANDERM4"]; # manual remove because of NaNs
  problem_list = filter_string_array(problem_list, remove_list);
  #test_problems(problem_list, 1)
  #CUTEstModel("GAUSS2")

  return problem_list
end

function get_problem_list(min_size::Int64, max_size::Int64)
  function filter_cutest(problem)
      min_size_ok = problem["constraints"]["number"] >= min_size && problem["variables"]["number"] >= min_size
      #max_size_ok = problem["constraints"]["number"] <= max_size && problem["variables"]["number"] <= max_size
      max_size_ok = problem["constraints"]["number"] + problem["variables"]["number"] <= max_size

      correct_size = min_size_ok && max_size_ok

      regular = problem["derivative_order"] >= 2 && problem["regular"] == true
      if correct_size && regular
          return true
      else
        return false
      end
  end

  return get_problem_list(filter_cutest)
end

function get_problem_list(min_ncon::Int64, max_ncon::Int64, min_nvar::Int64, max_nvar::Int64)
  function filter_cutest(problem)
      min_size_ok = problem["constraints"]["number"] >= min_ncon && problem["variables"]["number"] >= min_nvar
      max_size_ok = problem["constraints"]["number"] <= max_ncon && problem["variables"]["number"] <= max_nvar
      #max_size_ok = problem["constraints"]["number"] + problem["variables"]["number"] <= max_size

      correct_size = min_size_ok && max_size_ok

      regular = problem["derivative_order"] >= 2 && problem["regular"] == true
      if correct_size && regular
          @show problem["constraints"]["number"], problem["variables"]["number"]
          return true
      else
        return false
      end
  end

  return get_problem_list(filter_cutest)
end

function select_problem_with_sparse_rows(problem_list::Array{String,1},max_density::Int64)
    sparse_names = Array{String,1}()
    dense_names = Array{String,1}()
    for problem_name in problem_list
        nlp_raw = CUTEstModel(problem_name)

        J = jac(nlp_raw, nlp_raw.meta.x0)
        if densest_row(J) <= max_density
          push!(sparse_names, problem_name)
        else
          push!(dense_names, problem_name)
        end

        finalize(nlp_raw)
    end

    return sparse_names, dense_names
end

function default_list(sparsify=false)
    problem_list = get_problem_list(50,10000)
    # only run problems with max row density 1000.
    if sparsify
      sparse_names, dense_names = select_problem_with_sparse_rows(problem_list, 1000)
      problem_list = sparse_names
      println("these problems are included ...")
      @show sparse_names
      println("these problems are excluded ...")
      @show dense_names
    end

    return problem_list
end
