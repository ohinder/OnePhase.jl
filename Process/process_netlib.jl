using MAT, JuMP, Ipopt
include("../include.jl")

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

function run_infeas_netlib(test_name, my_par, solve_func!)
    summary = Dict{String, problem_summary}()
    problem_list = get_file_names("../netlib-infeas",".mat")

    if_mkdir("../results/$test_name")
    if_mkdir("../results/$test_name/log")
    if_mkdir("../results/$test_name/jld")

    for problem_name in problem_list
        if problem_name != "lpi_cplex2"
          m = read_lp2("../netlib-infeas", problem_name)
          nlp = MathProgNLPModel(m);
          solve_func!(nlp, problem_name, test_name, my_par, summary)
        end
    end
end

my_par = Class_parameters()
run_infeas_netlib("one_phase/infeas-lp-test", my_par, one_phase_solve!)
run_infeas_netlib("ipopt/infeas-lp-test", my_par, ipopt_solve!)
