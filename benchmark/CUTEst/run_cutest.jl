include("../benchmark.jl")
#using OnePhase, advanced_timer
using advanced_timer
using NLPModelsIpopt

if OnePhase.USE_HSL
	OnePhase.loadHSL("../../src/linear_system_solvers/")
end

function run_cutest_problems_on_solver(test_name::String, print_level::Int64, tol::Float64, max_iter::Int64, max_cpu_time::Float64, linear_solver::String, min_ncon::Int64, max_ncon::Int64, min_nvar::Int64, max_nvar::Int64)
	problems = get_problem_list(min_ncon, max_ncon, min_nvar, max_nvar)
	nlp_scaling_method = "none"
	bound_relax_factor =  0.0
	acceptable_iter = 999999
	run_cutest_problems_on_solver(problems, test_name, print_level, tol, max_iter, max_cpu_time, nlp_scaling_method, bound_relax_factor, acceptable_iter)
end

function run_cutest_problems_on_solver(problems::Array{String,1}, test_name::String, print_level::Int64, tol::Float64, max_iter::Int64, max_cpu_time::Float64, nlp_scaling_method::String, bound_relax_factor::Float64, acceptable_iter::Int64, linear_solver::String="")
    # if_mkdir("../results/$test_name")
    # if_mkdir("../results/$test_name/log")
	if_mkdir("$test_name")
    if_mkdir("$test_name/log")

	if isfile("$(test_name)/summary.jld")
	  # summary = load("../results/$(test_name)/summary.jld","summary")
	  summary = load("$(test_name)/summary.jld","summary")
	else
	  summary = Dict{String, problem_summary2}()
	end

	already_solved_problems = keys(summary)

    for problem_name in problems
		if problem_name in already_solved_problems
			 println("$problem_name already solved")
		else
			println("RUNNING $problem_name")
			ORG_STDOUT = stdout
	        # file = open("../results/$(test_name)/log/$(problem_name).txt", "w")
			file = open("$(test_name)/log/$(problem_name).txt", "w")
	        redirect_stdout(file)
	        summary[problem_name] = problem_summary2()
	        #tic()
			start_time = time()
	        try
				try
					nlp_raw = CUTEstModel(problem_name)
				catch(e)
					println("$problem_name failed to load")
					throw(e)
                end

				global t = 0
				function callback(alg_mod, iter_count, args...)
	           		global t = iter_count
	           		return true
	         	end
				stats = nothing
				if linear_solver != ""
					stats = ipopt(nlp_raw, print_level=print_level, tol=tol, max_iter=max_iter, max_cpu_time=max_cpu_time, nlp_scaling_method=nlp_scaling_method, bound_relax_factor=bound_relax_factor, acceptable_iter=acceptable_iter, linear_solver = linear_solver)
				else
					stats = ipopt(nlp_raw, print_level=print_level, tol=tol, max_iter=max_iter, max_cpu_time=max_cpu_time, nlp_scaling_method=nlp_scaling_method, bound_relax_factor=bound_relax_factor, acceptable_iter=acceptable_iter)
				end
				status = stats.status
	            summary[problem_name].status = status;
	            summary[problem_name].it_count = stats.iter;
		    	x = copy(stats.solution)
	            set_cutest_info_ipopt!(summary[problem_name], stats, nlp_raw, x)
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

          	# save("../results/$(test_name)/summary.jld","summary",summary)
		  	save("$(test_name)/summary.jld","summary",summary)

          	# summary_file = open("../results/$(test_name)/summary.txt", "w")
		  	summary_file = open("$(test_name)/summary.txt", "w")
          	write_summary(summary_file, summary)
          	close(summary_file)
    	end
	end
end


function run_cutest_problems_using_our_solver(test_name::String, par::OnePhase.Class_parameters, min_ncon::Int64, max_ncon::Int64, min_nvar::Int64, max_nvar::Int64)
	problems = get_problem_list(min_ncon, max_ncon, min_nvar, max_nvar)
	run_cutest_problems_using_our_solver(problems, test_name, par)
end

function run_cutest_problems_using_our_solver(problems::Array{String,1}, test_name::String, par::OnePhase.Class_parameters)

    # if_mkdir("../results/$test_name")
    # if_mkdir("../results/$test_name/log")
    # if_mkdir("../results/$test_name/jld")
	if_mkdir("$test_name")
    if_mkdir("$test_name/log")
    if_mkdir("$test_name/jld")

	# if isfile("../results/$(test_name)/summary.jld")
    if isfile("$(test_name)/summary.jld")
      # summary = load("../results/$(test_name)/summary.jld","summary")
	  summary = load("$(test_name)/summary.jld","summary")
    else
      summary = Dict{String, problem_summary2}()
    end

    already_solved_problems = keys(summary)

    # par_file = open("../results/$(test_name)/par.txt", "w")
	par_file = open("$(test_name)/par.txt", "w")
    OnePhase.write_pars(par_file, par)
    close(par_file)

    master_timer = class_advanced_timer()

    for problem_name in problems
      if problem_name in already_solved_problems
          println("$problem_name already solved")
      else
          println("RUNNING $problem_name")
          ORG_STDOUT = stdout
          # file = open("../results/$(test_name)/log/$(problem_name).txt", "w")
		  file = open("$(test_name)/log/$(problem_name).txt", "w")
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

              #tic()
	      start_time = time()

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
			  if OnePhase.USE_HSL
				  par.kkt.kkt_solver_type = :clever_symmetric
				  par.kkt.linear_solver_type = :HSL
			  end
              iter, status, history, t, err, timer = OnePhase.one_phase_solve(nlp_raw, par)

              master_timer = merge_timers(timer, master_timer)

              # save("../results/$(test_name)/jld/$(problem_name).jld","history",history, "timer", timer)
			  save("$(test_name)/jld/$(problem_name).jld","history",history, "timer", timer)

              set_info_me!(summary[problem_name], history, status, iter)
              #.it_count = t;
              summary[problem_name].number_variables = nlp_raw.meta.nvar
              summary[problem_name].number_constraints = nlp_raw.meta.ncon
              #summary[problem_name].total_time = toc();
              summary[problem_name].total_time = time() - start_time;
          catch(e)
              println("Uncaught error in algorithm!!!")
              @show e;

              if isa(e, OnePhase.Eval_NaN_error)
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

          # save("../results/$(test_name)/summary.jld","summary",summary, "pars", par)
		  save("$(test_name)/summary.jld","summary",summary, "pars", par)

          # summary_file = open("../results/$(test_name)/summary.txt", "w")
		  summary_file = open("$(test_name)/summary.txt", "w")
          write_summary(summary_file, summary)
          close(summary_file)

          # timer_file = open("../results/$(test_name)/timer.txt", "w")
		  timer_file = open("$(test_name)/timer.txt", "w")
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
    problem_list = get_problem_list(10,10000)
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
