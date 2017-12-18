include("../include.jl")

results = Dict{String, Dict{String,problem_summary}}()
results["one phase"] = load("../results/one_phase/infeas-test/summary.jld", "summary")
results["ipopt"] = convert_JuMP(load("../results/ipopt/infeas-test/summary.jld", "summary"))
overlapping_results = results;

##
## do exactly same thing as before ????
##
mode = :optimal

overlapping_results = all_status(overlapping_results,:primal_infeasible)
its, best, ratios, times = compute_its_etc(overlapping_results,MAX_IT=3000);

using PyPlot
plot_iteration_ratios(its, best, ratios)
