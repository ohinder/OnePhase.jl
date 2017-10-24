include("../include.jl")

results = Dict{String, Dict{String,problem_summary}}()
#results["IPOPT"] = load_results("other_solver_results/ipopt-results.txt");
#results["big_run"] = load("../results/OLD-august-21/big_run/summary.jld", "summary")
#results["IPOPT"] = convert_JuMP(load("../results/OLD-august-21/ipopt_test3/summary.jld", "summary"))

#results["IPOPT no perturb"] = convert_JuMP(load("../results/ipopt/large_no_perturb/summary.jld", "summary"))
#results["IPOPT"] = convert_JuMP(load("../results/ipopt/large/summary.jld", "summary"))
#results["static"] = load("../results/one_phase/large/summary.jld", "summary")
#results["large_dynamic"] = load("../results/one_phase/large_dynamic/summary.jld", "summary")
#results["dynamic"] = load("../results/one_phase/large_dynamic_1hr/summary.jld", "summary")
#results["one phase 2"] = load("../results/one_phase/plain/summary.jld", "summary")


#results["1"] = load("../results/one_phase/sept_1_corrections/summary.jld", "summary")
#results["2"] = load("../results/one_phase/sept_2_corrections/summary.jld", "summary")
#=
results["3"] = load("../results/one_phase/sept_3_corrections/summary.jld", "summary")
results["4"] = load("../results/one_phase/sept_4_corrections/summary.jld", "summary")
=#

# #=
#results["0"] = load("../results/one_phase/dual_ls_0/summary.jld", "summary")
#results["1"] = load("../results/one_phase/dual_ls_1/summary.jld", "summary")
#results["2"] = load("../results/one_phase/dual_ls_2/summary.jld", "summary")
# =#

# infeasible problems
#=
results["one phase"] = load("../results/one_phase/infeas-3/summary.jld", "summary")
results["IPOPT"] = convert_JuMP(load("../results/ipopt/infeas-3/summary.jld", "summary"))
=#

# #=
#results["max step"] = load("../results/one_phase/sept_max_step_stable/summary.jld", "summary")
#results["log barrier"] = load("../results/one_phase/sept_log_barrier_stable/summary.jld", "summary")
#results["filter"] = load("../results/one_phase/sept_3_corrections/summary.jld", "summary")
# =#

#=
results["regularizer"] = load("../results/one_phase/sept_no_regularizer/summary.jld", "summary")
results["no_regularizer"] = load("../results/one_phase/sept_3_corrections/summary.jld", "summary")
=#


#results["eta_affine"] = load("../results/one_phase/eta/affine/summary.jld", "summary")
#results["eta_mehrotra_stb"] = load("../results/one_phase/eta/mehrotra_stb/summary.jld", "summary")


#=
results["norm tol"] = load("../results/one_phase/norm_tol/summary.jld", "summary")
results["high tol"] = load("../results/one_phase/high_tol/summary.jld", "summary")
=#


#=
results["medium"] = convert_JuMP(load("../results/ipopt/medium/summary.jld", "summary"))
results["medium high tol"] = convert_JuMP(load("../results/ipopt/medium_high_tol/summary.jld", "summary"))
results["medium perturb"] = convert_JuMP(load("../results/ipopt/medium_perturb/summary.jld", "summary"))
results["medium perturb high tol"] = convert_JuMP(load("../results/ipopt/medium_perturb_high_tol/summary.jld", "summary"))
=#


# #=
#results["dynamic"] = load("../results/one_phase/sept_dynamic/summary.jld", "summary")
#results["static"] = load("../results/one_phase/sept_3_corrections/summary.jld", "summary")
# =#

# compare one phase and ipopt

#results["one phase"] = load("../results/one_phase/plain/summary.jld", "summary")
#results["ipopt"] = convert_JuMP(load("../results/ipopt/plain/summary.jld", "summary"))

#results["sept 2"] = load("../results/one_phase/sept_2_corrections/summary.jld", "summary")
#results["Oct21"] = load("../results/one_phase/Oct21/summary.jld", "summary")

#=
results["Oct21_penalty"] = load("../results/one_phase/Oct21_penalty/summary.jld", "summary")
results["Oct21_no_penalty"] = load("../results/one_phase/Oct21_no_penalty/summary.jld", "summary")
results["Oct21_penalty_full"] = load("../results/one_phase/Oct21_penalty_full/summary.jld", "summary")
results["Oct21_penalty_prox"] = load("../results/one_phase/Oct21_penalty_prox/summary.jld", "summary")
=#

#results["Oct22_a_penalty_only"] = load("../results/one_phase/Oct22/a_penalty_only/summary.jld", "summary")
#results["Oct22_a_penalty_only_no_prox"] = load("../results/one_phase/Oct22/a_penalty_only_no_prox/summary.jld", "summary")
results["Oct22_penalty"] = load("../results/one_phase/Oct22/penalty/summary.jld", "summary")
#results["Oct22_penalty_prox"] = load("../results/one_phase/Oct22/penalty_no_prox/summary.jld", "summary")




min_num = 0
#min_num = 100
max_num = 10000
function lrg_problems(problem)
    regular = problem["derivative_order"] >= 2 && problem["regular"] == true && problem["constraints"]["number"] >= 1
    min_size_ok = problem["variables"]["number"] >= min_num && problem["constraints"]["number"] >= min_num
    #max_size = problem["variables"]["number"] <= max_num && problem["constraints"]["number"] <= max_num
    max_size_ok = problem["constraints"]["number"] + problem["variables"]["number"] <= max_num

    if min_size_ok && max_size_ok && regular
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
#error_free_results = remove_errors(results, [:NaN_ERR,:Error])
#error_free_results = remove_errors(results, [:MAX_TIME])
#error_free_results = remove_errors(results, [:NaN_ERR, :MAX_TIME])
error_free_results = remove_errors(results, [:MAX_TIME])
overlapping_results = overlap(error_free_results)
#overlapping_results = restrict_to_set(overlapping_results,[:primal_infeasible])


####overlapping_results = restrict_to_set(overlapping_results,[:optimal,:primal_infeasible,:dual_infeasible])

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

println("kkt")
for (method_name, sum_data) in overlapping_results
    println(method_name, "=" ,tot(sum_data,[:optimal]), " out of ", length(sum_data))
end

println("unbounded")
for (method_name, sum_data) in overlapping_results
    println(method_name, "=" ,tot(sum_data,[:dual_infeasible]), " out of ", length(sum_data))
end

println("primal infeasible")
for (method_name, sum_data) in overlapping_results
    println(method_name, "=" ,tot(sum_data,[:primal_infeasible]), " out of ", length(sum_data))
end

println("failures")
for (method_name, sum_data) in overlapping_results
    println(method_name, "=" ,tot_failures(sum_data), " out of ", length(sum_data))
end

@show shared_failures(overlapping_results)
@show list_combined_success(overlapping_results)

#outcomes_table(overlapping_results)


println("failures")
for (method_name, sum_data) in overlapping_results
    println(method_name)#
    @show failure_causes(sum_data)
    #, "=" ,tot_failures(sum_data), " out of ", length(sum_data))
end



print_failure_problems(overlapping_results)
its, best, ratios, times = compute_its_etc(overlapping_results);


##
## PLOT ITERATION RATIO CURVES
##

using PyPlot

plot_iteration_ratios(its, best, ratios)
plot_iterations(its, best, ratios, 3000)



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
      both_kkt = true
      fbest = 1e14
      for method_name in method_list
        info = overlapping_results[method_name][problem_name]
        fval = info.fval
        con_vio = info.con_vio
        print(rd(con_vio))
        print(rd(fval))

        if info.status == :optimal
          if con_vio > 1e-6
            both_kkt = false
            println("CON VIOLATION!!")
          else
            @assert(con_vio < 1e-6)
            ftol = 1e-1
            ftol_scaled = ftol * (1.0 + abs(fbest))
            if fval <= fbest - ftol_scaled
              fbest = fval
              best = method_name
            elseif fbest - ftol_scaled < fval && fval <= fbest + ftol_scaled
              fbest = fval
              best = "tie"
            end
          end
        else
          both_kkt = false
        end
      end
      if both_kkt
        winners[best] += 1
      end
      print(pd(best))
      print("\n")
  end
end

winners

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
method_name2 = "IPOPT no perturb"
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
ax[:set_xlim]([1.0,2^8.0])
ax[:set_ylim]([min_y,1.0])

#ax[:xaxis][:ticker] = 0.5
#title("Comparsion on 45 CUTEst problems")
xlabel("times ratio to best solver")
ylabel("proportion of problems")

legend()
