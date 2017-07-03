include("include.jl")

results = Dict{String, Dict{String,problem_summary}}()
#results["IPOPT"] = load_results("other_solver_results/ipopt-results.txt");
#results["KNITRO"] = load_results("other_solver_results/knitro-results.txt");
#results["ME1"] = load("results/test6/summary.jld", "summary")
#results["ME1"] = load("results/test7/summary.jld", "summary")
#results["Mehotra2"] = load("results/mehotra_intial_point2/summary.jld", "summary")

#results["One Phase"] = load("results/mehotra_intial_point3/summary.jld", "summary")
#results["One Phase 2"] = load("results/mehotra_intial_point4/summary.jld", "summary")
#results["One Phase 3"] = load("results/mehotra_intial_point5/summary.jld", "summary")
#results["One Phase 4"] = load("results/mehotra_intial_point6-copy/summary.jld", "summary")
#results["One Phase 5"] = load("results/mehotra_intial_point7/summary.jld", "summary")

#results["One Phase 3"] = load("results/pars3/mu-fix/summary.jld", "summary")
#results["One Phase 4"] = load("results/pars3/mu-test1/summary.jld", "summary")
#results["One Phase 5"] = load("results/pars3/mu-test3/summary.jld", "summary")
#results["One Phase 6"] = load("results/pars3/mu-awesome/summary.jld", "summary")


#results["no inertia test"] = load("results/inertia_test/true/summary.jld", "summary")
#results["inertia test"] = load("results/inertia_test/false/summary.jld", "summary")


#results["One Phase 5"] = load("results/pars4/test5/summary.jld", "summary")
#results["One Phase 6"] = load("results/pars4/test6/summary.jld", "summary")
#results["One Phase 7"] = load("results/pars4/test7/summary.jld", "summary")
#results["One Phase 8"] = load("results/pars4/test8/summary.jld", "summary")
#results["One Phase 9"] = load("results/pars4/test9/summary.jld", "summary")
#results["One Phase 10"] = load("results/pars4/test10/summary.jld", "summary")
#results["One Phase 11"] = load("results/pars4/test11/summary.jld", "summary")
#results["One Phase 12"] = load("results/pars4/test12/summary.jld", "summary")
results["IPOPT"] = convert_JuMP(load("results/ipopt_test3/summary.jld", "summary"))
#results["One Phase none"] = load("results/pars4/none/summary.jld", "summary")
#results["One Phase true"] = load("results/pd-split/true/summary.jld", "summary")
#results["One Phase false"] = load("results/pd-split/false/summary.jld", "summary")
#results["One Phase ls3"] = load("results/ls/ls_true3/summary.jld", "summary")
#results["One Phase ls2"] = load("results/ls/ls_true2/summary.jld", "summary")
#results["One Phase ls1"] = load("results/ls/ls_true/summary.jld", "summary")
#results["One Phase ls false"] = load("results/ls/ls_false/summary.jld", "summary")
#results["new approach"] = load("results/new_approach/summary.jld", "summary")
#results["new approach 2"] = load("results/new_approach2/summary.jld", "summary")
#results["new approach 3"] = load("results/new_approach3/summary.jld", "summary")
#results["new approach no prox"] = load("results/new_approach_no_prox/summary.jld", "summary")
#results["new approach lrg dual"] = load("results/new_approach_no_prox/summary.jld", "summary")
#results["new approach"] = load("results/new_approach_latest/summary.jld", "summary")
#results["new approach 2"] = load("results/new_approach_latest/summary.jld", "summary")
#results["new approach 2"] = load("results/new_approach_latest_working/summary.jld", "summary")
#results["new approach 3"] = load("results/new_approach_latest2/summary.jld", "summary")

#results["new approach 4"] = load("results/new_approach_latest3/summary.jld", "summary")
results["big_run"] = load("results/big_run/summary.jld", "summary")

function restrict_problems(results, problem_list)
  new_results = Dict{String,Dict{String,problem_summary}}()
  for (method_name, sum_data) in results
      new_results[method_name] = Dict()
      for (problem_name, info) in sum_data
        if problem_name in problem_list
          new_results[method_name][problem_name] = info
        end
      end
  end
  return new_results
end

function lrg_problems(problem)
    regular = problem["derivative_order"] >= 2 && problem["regular"] == true
    correct_size = problem["variables"]["number"] + problem["constraints"]["number"] >= 500 #&& problem["variables"]["number"] + problem["constraints"]["number"] <= 1600
    if correct_size && regular
        return true
    else
      return false
    end
end

function real_problems(problem)
    if problem["origin"] == "real" || problem["origin"] == "modelling"
        return true
    else
        return false
    end
end

problem_list = CUTEst.select(custom_filter=lrg_problems);
problem_list = convert(Array{String,1},problem_list);

results = restrict_problems(results, problem_list)

#results["stable_first/false"] = load("results/stable_first/false/summary.jld", "summary")
#results["stable_first/true"] = load("results/stable_first/true/summary.jld", "summary")

if true
error_free_results = remove_errors(results, [:NaN_ERR]) #, :MAX_TIME])#, :ERR])
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
its, best, ratios, times = compute_its_etc(overlapping_results)


##
## PLOT ITERATION RATIO CURVES
##

using PyPlot

min_y = 1.0

y_vals = collect(1:length(best)) / length(best)
for (method_name, val) in its
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
ax[:set_xlim]([10.0,my_par.max_it])
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

winners = Dict()
winners["tie"] = 0
winners["fail"] = 0
begin
  for method_name in method_list
      print(pd(""), pd(method_name))
      winners[method_name] = 0
  end
  print("\n")

  for problem_name in problem_list
      print(pd(problem_name))
      best = "fail"
      fbest = 1e14
      for method_name in method_list
        fval = overlapping_results[method_name][problem_name].fval
        con_vio = overlapping_results[method_name][problem_name].con_vio
        print(rd(con_vio))
        print(rd(fval))
        if overlapping_results[method_name][problem_name].con_vio < 1e-5
          ftol = 1e-2
          ftol_scaled = ftol * (1.0 + abs(fbest))
          if fval <= fbest - ftol_scaled
            fbest = fval
            best = method_name
          elseif fbest - ftol_scaled < fval && fval <= fbest + ftol_scaled
            fbest = fval
            best = "tie"
          end
        end
      end
      winners[best] += 1
      print(pd(best))
      print("\n")
  end
end

function best_fvals(method_name1, method_name2)
    better_fval = Dict()
    succeed = Dict()

    feas_tol = 1e-4
    fval_tol = 1e-1

    better_fval[method_name1] = 0
    better_fval[method_name2] = 0

    succeed[method_name1] = 0
    succeed[method_name2] = 0

    for problem_name in problem_list
      fval1 = overlapping_results[method_name1][problem_name].fval
      con_vio1 = overlapping_results[method_name1][problem_name].con_vio

      fval2 = overlapping_results[method_name2][problem_name].fval
      con_vio2 = overlapping_results[method_name2][problem_name].con_vio

      if con_vio1 < feas_tol && con_vio2 < feas_tol
        fbest = min(fval1,fval2)
        ftol_scaled = fval_tol * (1.0 + abs(fbest))

        if fval1 < fval2 - ftol_scaled
          better_fval[method_name1] += 1
        elseif fval2 < fval1 - ftol_scaled
          better_fval[method_name2] += 1
        end
      end

      if con_vio1 < feas_tol
        succeed[method_name1] += 1
      end

      if con_vio2 < feas_tol
        succeed[method_name2] += 1
      end
    end

    return better_fval, succeed
end

method_name1 = "IPOPT"
method_name2 = "new approach 4"
best_fvals(method_name1, method_name2)



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
ax[:set_xlim]([1.0,2^5.0])
ax[:set_ylim]([min_y,1.0])

#ax[:xaxis][:ticker] = 0.5
title("Comparsion on 119 CUTEst problems")
#title("Comparsion on 45 CUTEst problems")
xlabel("iteration ratio to best solver")
ylabel("proportion of problems")

legend()
