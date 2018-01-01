include("../benchmark.jl")
include("create_report.jl")

overlapping_results = get_CUTEst_results()
#dual_inf_free_results = remove_errors(overlapping_results, [:dual_infeasible,:primal_infeasible])
#overlapping_opt_results = overlap(dual_inf_free_results)
#overlapping_opt_results = restrict_to_set(overlapping_opt_results,[:optimal])
#overlapping_opt_results = overlap(overlapping_opt_results)

####
#### Lets do an overall summary of the results
####

# i.e. create Table 3 through 5.

# compute median number of iterations


##### Comparison on time to find a KKT point
#####
##### Lets take a look at the problems for which
##### - both solvers return a KKT point with the same objective
#####
#####

using PyPlot
PyPlot.close()
its, best, ratios, times = compute_its_etc(overlapping_results,MAX_IT=3000);
plot_iterations(its, best, ratios, 3000, line_styles, line_colors)
savefig("$folder/iter_curve.pdf")
PyPlot.close()

plot_iteration_ratios(its, best, ratios, line_styles, line_colors)
savefig("$folder/iter_ratios.pdf")
PyPlot.close()


res = all_status(overlapping_results, mode)
its, best, ratios, times = compute_its_etc(res,MAX_IT=3000);

println("number of overlapping results = ", length(res["ipopt"]))

plot_iteration_ratios(its, best, ratios)
savefig("$folder/opt_iter_ratios.pdf")
PyPlot.close()

##
## PLOT TIME CURVES
##
best = best_its(times)

min_y = 1.0

ratios = Dict()
y_vals = collect(1:length(best)) / length(best)
for (method_name, val) in its
  ratios[method_name] = times[method_name] ./ best
  semilogx(sort(ratios[method_name]), y_vals, label=method_name, basex=2)
  min_y = min(min_y,sum(ratios[method_name] .== 1.0) / length(best))
end
ax = gca()
ax[:set_xlim]([1.0,2^8.0])
ax[:set_ylim]([min_y,1.0])

#ax[:xaxis][:ticker] = 0.5
#title("Comparsion on 45 CUTEst problems")
xlabel("times ratio to best solver")
ylabel("proportion of problems")

legend()
