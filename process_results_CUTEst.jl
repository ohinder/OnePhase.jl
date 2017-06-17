include("include.jl")

results = Dict{String, Dict{String,problem_summary}}()
#results["IPOPT"] = load_results("other_solver_results/ipopt-results.txt");
#results["KNITRO"] = load_results("other_solver_results/knitro-results.txt");
#results["ME1"] = load("results/test6/summary.jld", "summary")
#results["ME1"] = load("results/test7/summary.jld", "summary")
#results["Mehotra2"] = load("results/mehotra_intial_point2/summary.jld", "summary")

#results["One Phase"] = load("results/mehotra_intial_point3/summary.jld", "summary")
#results["One Phase 2"] = load("results/mehotra_intial_point4/summary.jld", "summary")
results["One Phase 3"] = load("results/mehotra_intial_point5/summary.jld", "summary")
results["One Phase 4"] = load("results/mehotra_intial_point6/summary.jld", "summary")

results["IPOPT"] = convert_JuMP(load("results/ipopt_test2/summary.jld", "summary"))

#results["grad scaled mu"] = load("results/inside2/summary.jld", "summary")

#results["ME3"] = load("results/par1/mu_ratio-0.01/summary.jld", "summary")
#results["ME4"] = load("results/par1/mu_ratio-1.0/summary.jld", "summary")
#results["ME5"] = load("results/par1/mu_ratio-100.0/summary.jld", "summary")
#results["mehotra"] = load("results/par2/mehotra-no-satisfy/summary.jld", "summary")
#results["mehotra"] = load("results/par2/mehotra-no-satisfy/summary.jld", "summary")
#results["mu 0.01"] = load("results/par_hess1/mu_ratio-0.01/summary.jld", "summary")
#results["mu 1.0"] = load("results/par_hess1/mu_ratio-1.0/summary.jld", "summary")
#results["mu 100.0"] = load("results/par_hess1/mu_ratio-100.0/summary.jld", "summary")

if false
error_free_results = remove_errors(results, [:NaN_ERR, :ERR])
overlapping_results = overlap(error_free_results)
elseif false
overlapping_results = overlap(results)
infeas_results = select(overlapping_results, [:primal_infeasible])
overlapping_results = union_results(infeas_results, results)
show_results(overlapping_results)
else
overlapping_results = overlap(results)
end

##
## FAILURES
##

for (method_name, sum_data) in overlapping_results
    println(method_name, "=" ,tot_failures(sum_data), " out of ", length(sum_data))
end

list_failures(overlapping_results)

its = Dict{String,Array{Int64,1}}()
for (method_name, sum_data) in overlapping_results
    its[method_name] = iteration_list(sum_data);
end



##
## PLOT ITERATION CURVES
##
best = best_its(its)

using PyPlot

min_y = 1.0

ratios = Dict()
y_vals = collect(1:length(best)) / length(best)
for (method_name, val) in its
  ratios[method_name] = its[method_name] ./ best
  semilogx(sort(ratios[method_name]), y_vals, label=method_name, basex=2)
  min_y = min(min_y,sum(ratios[method_name] .== 1.0) / length(best))
end
ax = gca()
ax[:set_xlim]([1.0,2^5.0])
ax[:set_ylim]([min_y,1.0])

#ax[:xaxis][:ticker] = 0.5
title("Comparsion on 119 CUTEst problems")
#title("Comparsion on 45 CUTEst problems")
xlabel("iteration ratio to best solver")
ylabel("proportion of problems")

legend()


min_y = 1.0

y_vals = collect(1:length(best)) / length(best)
for (method_name, val) in its
  semilogx(sort(its[method_name]), y_vals, label=method_name, basex=10)
  min_y = min(min_y,sum(its[method_name] .== 1.0) / length(best))
end
ax = gca()
ax[:set_xlim]([10.0,10^4.0])
ax[:set_ylim]([min_y,1.0])

#ax[:xaxis][:ticker] = 0.5
title("Comparsion on 119 CUTEst problems")
#title("Comparsion on 45 CUTEst problems")
xlabel("iterations")
ylabel("proportion of problems")

legend()


##
## FUNCTION VALUES
##
problem_list = collect(keys(first(overlapping_results)[2]))
method_list = collect(keys(overlapping_results))

for method_name in method_list
    print(pd(method_name))
end
print("\n")

for problem_name in problem_list
    print(pd(problem_name))
    for method_name in method_list
      print(rd(overlapping_results[method_name][problem_name].con_vio))
      print(rd(overlapping_results[method_name][problem_name].fval))
    end
    print("\n")
end

##
## DUAL sequences
##

#=
problem_list = collect(keys(first(overlapping_results)[2]))
method_list = keys(overlapping_results)

dual_maxes = Dict{String,Array{Float64,1}}()
for folder_name in method_list
  dual_maxes[folder_name] = Array{Float64,1}()
  for problem_name in problem_list
      results = load("results/$folder_name/jld/$problem_name.jld", "history")
      dual_hist = compute_dual_hist(results)
      push!(dual_maxes[folder_name], maximum(dual_hist))
  end
end

best_duals = best_its(dual_maxes)

min_y = 1.0

ratios = Dict()
y_vals = collect(1:length(best_duals)) / length(best_duals)
for (method_name, val) in dual_maxes
  #ratios[method_name] = dual_maxes[method_name] ./ best_duals
  semilogx(sort(dual_maxes[method_name]), y_vals, label=method_name, basex=10)
  #min_y = min(min_y,sum(ratios[method_name] .== 1.0) / length(best_duals))
end

xlabel("maximum dual value")
ylabel("proportion of problems")

legend()=#
