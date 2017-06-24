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
          summary[problem_name] = problem_summary()
          start_time = time()

          nlp_raw = false

          try
              try
                nlp_raw = CUTEstModel(problem_name)
              catch(e)
                println("$problem_name failed to load")
                throw(e)
              end

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
  regular = problem["derivative_order"] >= 2 && problem["regular"] == true
    # large
    #correct_size = 50 <= problem["variables"]["number"] + problem["constraints"]["number"] && problem["constraints"]["number"] >= 10 && problem["variables"]["number"] + problem["constraints"]["number"] <= 2000
    # medium lrg
    #10 <= problem["variables"]["number"] + problem["constraints"]["number"] &&
    correct_size = problem["constraints"]["number"] >= 1 && problem["variables"]["number"] + problem["constraints"]["number"] <= 2000
    # medium
    #correct_size = 50 <= problem["variables"]["number"] + problem["constraints"]["number"] && problem["constraints"]["number"] >= 10 && problem["variables"]["number"] <= 1000 && problem["constraints"]["number"] <= 1000
    # small
    #correct_size = 100 <= problem["variables"]["number"] + problem["constraints"]["number"] && problem["constraints"]["number"] >= 10 && problem["variables"]["number"] + problem["constraints"]["number"] <= 300
    if correct_size && regular
        return true
    else
      return false
    end
end

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

#problem_list = CUTEst.select(max_var=100, max_var=1000, min_con=100, max_con=3000)
problem_list = CUTEst.select(custom_filter=filter_cutest);
problem_list = convert(Array{String,1},problem_list);
remove_list = ["MISRA1D","OSBORNE1","LANCZOS2","MEYER3NE","ROSZMAN1","INTEQNE","MGH10",
"ECKERLE4","RAT43","MISRA1A","KOWOSBNE","MOREBVNE",
"JENSMPNE","MISRA1C","VARDIMNE","SANTA","BARDNE","MGH10S",
"THURBER","MGH17S","BOX3NE","PENLT1NE","GAUSS2","GULFNE","BA-L1SP","ENSO","CHWIRUT1", "CHWIRUT2","NELSON","HAHN1","KIRBY2","MUONSINE","OSBORNE2","BENNETT5","BA-L1","GAUSS3",
"GAUSS1","LANCZOS3","BIGGS6NE","LANCZOS1","VANDANIUMS","MISRA1B","MGH09","MGH17","WATSONNE",
"DMN37142LS","VESUVIOU","PENLT2NE", "DMN37143LS","INTEQNELS","BROYDN3DLS","DMN15332LS",
"DMN15102LS","BROYDNBDLS","VESUVIA","BA-L1SPLS","BA-L1SPLS","SANTALS","VESUVIO","DMN15333LS","ARGTRIGLS","RAT42","SSINE","LSC1","LSC2","BOXBOD","POWELLSE","FREURONE","DANWOOD","HELIXNE"];
problem_list = filter_string_array(problem_list, remove_list);
test_problems(problem_list,1)
#CUTEstModel("GAUSS2")

if false
problem_list = ["PT"]
folder_name = "test_run"
if_mkdir("results/$folder_name")
run_cutest_problems_using_our_solver(problem_list, folder_name, my_par)
end

if false
    folder_name = "stable_first"
    if_mkdir("results/$folder_name")

    subfolder_name = "$folder_name/true"
    if_mkdir("results/$subfolder_name")
    my_par.stb_before_agg = true
    my_par.max_it_corrections = 2
    run_cutest_problems_using_our_solver(problem_list, subfolder_name, my_par)

    subfolder_name = "$folder_name/false"
    if_mkdir("results/$subfolder_name")
    my_par.stb_before_agg = false
    my_par.max_it_corrections = 3
    run_cutest_problems_using_our_solver(problem_list, subfolder_name, my_par)
end

if false
    main_folder_name = "pars4"
    if_mkdir("results/$main_folder_name")

    style_list = [:none, :test5, :test6, :test7, :test8, :test9, :test10, :test11, :test12]
    for mu_style in style_list
      path = "$main_folder_name/$mu_style"
      if_mkdir("results/$path")
      my_par.adaptive_mu = mu_style
      run_cutest_problems_using_our_solver(problem_list, path, my_par)
    end

    #subfolder_name = "$folder_name/mu-fix"
    #if_mkdir("results/$subfolder_name")
    #my_par.adaptive_mu = :none
    #run_cutest_problems_using_our_solver(problem_list, subfolder_name, my_par)

    #folder_name = "par2/mehotra-no-satisfy"
    #if_mkdir("results/$folder_name")
    #my_par.start_satisfying_bounds = false
    #run_cutest_problems_using_our_solver(problem_list, "$folder_name", my_par)
end

if true
  main_folder_name = "ls"
  if_mkdir("results/$main_folder_name")

  #for mu_style in [true; false]
  path = "$main_folder_name/ls_true3"
  if_mkdir("results/$path")
  my_par.move_primal_seperate_to_dual = true #mu_style
  my_par.dual_ls = true
  run_cutest_problems_using_our_solver(problem_list, path, my_par)
  #end
end



if false
    main_folder_name = "inertia_test"
    if_mkdir("results/$main_folder_name")


    for mu_style in [true; false]
      path = "$main_folder_name/$mu_style"
      if_mkdir("results/$path")
      my_par.inertia_test = mu_style
      run_cutest_problems_using_our_solver(problem_list, path, my_par)
    end

    #subfolder_name = "$folder_name/mu-fix"
    #if_mkdir("results/$subfolder_name")
    #my_par.adaptive_mu = :none
    #run_cutest_problems_using_our_solver(problem_list, subfolder_name, my_par)

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
