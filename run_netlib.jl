include("include.jl")
include("read_lp.jl")

using Ipopt


function run_netlib_problems_using_IPOPT(problem_list::Array{String,1}, test_name::String, bound_relax_factor::Float64)
    max_iter = 100

    summary = Dict{String, problem_summary}()

    folder_name = "results/$test_name"
    jld_folder = "$folder_name/jld"
    log_folder = "$folder_name/log"
    if_mkdir(folder_name)
    if_mkdir(jld_folder)
    if_mkdir(log_folder)

    for problem_name in problem_list
        start_time = time()
        println("RUNNING ", problem_name)

        ORG_STDOUT = STDOUT
        file = open("$log_folder/$(problem_name).txt", "w")
        #redirect_stdout(file)

        A, b, c, lbounds, ubounds = read_lp(problem_name)

        solver = IpoptSolver(print_level=5, max_iter=max_iter)
        m = read_lp_into_JuMP(problem_name)

        setsolver(m, solver)
        JuMP.build(m)
        status = solve(m)
        @show status


        if status == :Optimal
          hist = run_IPOPT(A, b, c, lbounds, ubounds, max_iter, bound_relax_factor)
          #run_IPOPT(m, max_iter, bound_relax_factor)
          save("$(jld_folder)/$(problem_name).jld","history",hist)

          summary[problem_name] = problem_summary()
          summary[problem_name].status = status;
          summary[problem_name].total_time = time() - start_time
          set_info_me!(summary[problem_name], hist)

          save("$(folder_name)/summary.jld","summary",summary)

          summary_file = open("$(folder_name)/summary.txt", "w")
          write_summary(summary_file, summary)
          close(summary_file)
        end

        redirect_stdout(STDOUT)
        close(file)
    end
end

function run_netlib_problems_using_our_solver(problem_list::Array{String,1}, test_name::String, par::Class_parameters)
    summary = Dict{String, problem_summary}()

    folder_name = "results/$test_name"
    jld_folder = "$folder_name/jld"
    log_folder = "$folder_name/log"
    if_mkdir(folder_name)
    if_mkdir(jld_folder)
    if_mkdir(log_folder)

    for problem_name in problem_list
        start_time = time()
        println("RUNNING ", problem_name)
        jump_model = read_lp_into_JuMP(problem_name, 0.0, true) #read_lp_into_Class_QP(problem_name);
        nlp_raw = MathProgNLPModel(jump_model)

        ORG_STDOUT = STDOUT
        file = open("$log_folder/$(problem_name).txt", "w")
        redirect_stdout(file)

        nlp = Class_CUTEst(nlp_raw)
        timer = class_advanced_timer()
        start_advanced_timer(timer)
        start_advanced_timer(timer, "INIT")
        intial_it = init(nlp, my_par, timer)
        pause_advanced_timer(timer, "INIT")
        @assert(is_feasible(intial_it, my_par.comp_feas))
        iter, status, hist, t, err = one_phase_IPM(intial_it, my_par, timer);
        pause_advanced_timer(timer)

        print_timer_stats(timer)

        finalize(nlp_raw)

        save("$(jld_folder)/$(problem_name).jld","history",hist)
        summary[problem_name] = problem_summary()
        summary[problem_name].status = status
        set_info_me!(summary[problem_name], hist)
        summary[problem_name].total_time = time() - start_time

        redirect_stdout(ORG_STDOUT)
        close(file)

        save("$(folder_name)/summary.jld", "summary", summary, "pars", par)

        summary_file = open("$(folder_name)/summary.txt", "w")
        write_summary(summary_file, summary)
        close(summary_file)
    end
end

###
### LOAD FILES
###


problem_list = load_netlib(Int(1e4))

function pars_well_behaved!(pars::Class_parameters)
    pars.adaptive_mu = :none
    pars.max_it = 100
    pars.aggressive_reduct_factors = Reduct_affine()
end

function pars_fast_primal!(pars::Class_parameters)
    pars.adaptive_mu = :none
    pars.max_it = 100
    pars.aggressive_reduct_factors = Class_reduction_factors(0.0,0.0,0.2)
end

begin
  if_mkdir("results/netlib")

  if true
  test_name = "netlib/fast-primal"
  pars_fast_primal!(my_par);
  run_netlib_problems_using_our_solver(problem_list, test_name, my_par)
  end

  if true
  test_name = "netlib/well-behaved"
  pars_well_behaved!(my_par);
  run_netlib_problems_using_our_solver(problem_list, test_name, my_par)
  end


  if true
  test_name = "netlib/ipopt-perturb"
  run_netlib_problems_using_IPOPT(problem_list, test_name, 1e-8)
  end

  if true
  # without increasing the size of the constraints
  test_name = "netlib/ipopt-no-perturb"
  run_netlib_problems_using_IPOPT(problem_list, test_name, 0.0)
  end
end

if false
include("include.jl")
jump_model = read_lp_into_JuMP("PEROLD", 0.0, true);
nlp_raw = MathProgNLPModel(jump_model);
nlp = Class_CUTEst(nlp_raw);
timer = class_advanced_timer();

start_advanced_timer(timer);
start_advanced_timer(timer, "INIT");
intial_it = init(nlp, my_par, timer);
pause_advanced_timer(timer, "INIT");
@assert(is_feasible(intial_it, my_par.comp_feas));
iter, status, hist, t, err = one_phase_IPM(intial_it, my_par, timer);
pause_advanced_timer(timer);

print_timer_stats(timer);
end
