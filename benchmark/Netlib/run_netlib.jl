using MAT, JuMP, Ipopt
include("../benchmark.jl")

function read_lp2(dir, name)
    lp_data = matopen("$(dir)/$(name).mat")

    problem = read(lp_data,"Problem")
    A = problem["A"]
    b = problem["b"][:]
    c = problem["aux"]["c"][:]
    lb = problem["aux"]["lo"][:]
    ub = problem["aux"]["hi"][:]
    n = length(c)

    m = Model()
    @variable(m, lb[i] <= x[i=1:n] <= ub[i])
    @NLobjective(m, Min, sum(c[i] * x[i] for i = 1:n))
    @constraint(m, A * x .== b)


    return m
end

function get_file_names(dir::String,ext::String)
  name_list = readdir(dir)
  cleaned_name_list = []
  for name_full in name_list
    try
      if name_full[(end-3):end] == ext
          name = name_full[1:(end-4)]
      end
      cleaned_name_list = [cleaned_name_list; name]
      catch(e)
        @show e
        println("ERROR in running " * name_full)
      end
  end

  return cleaned_name_list
end

function run_infeas_netlib(test_name, my_par, solve_func!)
    summary = Dict{String, problem_summary2}()
    dir = "../data/netlib-infeas"
    problem_list = get_file_names(dir,".mat")

    if_mkdir("../results/$test_name")
    if_mkdir("../results/$test_name/log")
    if_mkdir("../results/$test_name/jld")

    for problem_name in problem_list
        if problem_name != "lpi_cplex2"
          m = read_lp2(dir, problem_name)
          nlp = MathProgNLPModel(m);
          #if nlp.meta.ncon < 2000
          summary[problem_name] = solve_func!(nlp, problem_name, test_name, my_par)
          #end
        end
        write_summary("../results/$(test_name)/summary.txt", summary)
        save("../results/$(test_name)/summary.jld","summary",summary, "pars", my_par)
    end
end

my_par = Class_parameters()
run_infeas_netlib("one_phase/Jan4/infeas-netlib", my_par, one_phase_run_and_store)
run_infeas_netlib("ipopt/Jan4/infeas-netlib", my_par, ipopt_solve_run_and_store)
