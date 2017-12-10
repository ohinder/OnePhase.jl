include("../include.jl")
include("create_report.jl")
overlapping_results = get_CUTEst_results()
dual_inf_free_results = remove_errors(overlapping_results, [:dual_infeasible,:primal_infeasible])
overlapping_opt_results = overlap(dual_inf_free_results)
overlapping_opt_results = restrict_to_set(overlapping_opt_results,[:optimal])
overlapping_opt_results = overlap(overlapping_opt_results)

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

function create_opt_res(res::Dict{String, Dict{String,problem_summary}}, method_list; f_TOL = 1e-1)
    opt_res = Dict{String, Dict{String,problem_summary}}()
    for method in method_list
        opt_res[method] = Dict{String,problem_summary}()
    end

    problem_list = keys(first(res)[2])
    for problem_name in problem_list
      same_fval = true
      first = true
      fval_others = Inf

      for method in method_list
         fval = res[method][problem_name].fval
         if res[method][problem_name].status == :optimal
           if first || abs(fval_others - fval) / (1.0 + max(abs(fval),abs(fval_others))) < f_TOL
             fval_others = min(fval_others, fval)
             first = false
           else
             same_fval = false
           end
         else
           same_fval = false
         end
      end

      if same_fval #|| some_method_failed
          for method in method_list
             opt_res[method][problem_name] = res[method][problem_name]
          end
      end
    end
    return opt_res
end

method_list = keys(overlapping_opt_results)
opt_res = create_opt_res(overlapping_opt_results, method_list, f_TOL = 1e-2)

its, best, ratios, times = compute_its_etc(opt_res,MAX_IT=3000);

using PyPlot

plot_iteration_ratios(its, best, ratios)
savefig("output/opt_iter_ratios.pdf")

plot_iterations(its, best, ratios, 3000)
savefig("output/opt_iter_curve.pdf")

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
