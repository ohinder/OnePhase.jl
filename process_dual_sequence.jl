include("include.jl")

summary_netlib = Dict{String, Dict{String,problem_summary}}()
summary_netlib["netlib"] = load("results/netlib/summary.jld", "summary")
summary_netlib["netlib-ipopt"] = convert_JuMP(load("results/netlib-ipopt/summary.jld", "summary"))


for (method_name, sum_data) in summary_netlib
    println(method_name, "=" ,tot_failures(sum_data), " out of ", length(sum_data))
end


#problem_name = "BANDM"
#hist =

function compute_dual_hist(results)
    dual_hist = Array{Float64,1}()
    t = 0
    for it_hist in results
      if it_hist.t > t
        push!( dual_hist, it_hist.y_norm )
        t += 1
      end
    end
    #dual_maxes[folder_name][problem_name] = maximum(dual_hist)

    return dual_hist
end

function get_lists(summary)
    problem_list = collect(keys(first(summary)[2]))
    method_list = keys(summary)

    return problem_list, method_list
end

problem_list, method_list = get_lists(summary_netlib)
#problem_list = collect(keys(summary_netlib["me"]))
#method_list = ["netlib","netlib-ipopt"]

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

using PyPlot

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

legend()

# IPOPT 11/29 failures
# us 4/29 failures




con_vio = Dict{String,Array{Float64,1}}()
for folder_name in method_list
  con_vio[folder_name] = Array{Float64,1}()
  for problem_name in problem_list
      results = load("results/$folder_name/jld/$problem_name.jld", "history")
      push!(con_vio[folder_name], results[end].primal_residual)
  end
end

ratios = Dict()
y_vals = collect(1:length(con_vio["netlib"])) / length(con_vio["netlib"])
for (method_name, val) in dual_maxes
  #ratios[method_name] = dual_maxes[method_name] ./ best_duals
  semilogx(sort(con_vio[method_name]), y_vals, label=method_name, basex=10)
  #min_y = min(min_y,sum(ratios[method_name] .== 1.0) / length(best_duals))
end

xlabel("con violation")
ylabel("proportion of problems")

legend()


#=
for problem_name in problem_list
    for folder_name in folder_list
      dual_maxes[]
    end
end=#




hist = Dict{String, Dict{String,abstract_alg_history}}()
hist["me"] = dual_hist["me"]
hist["IPOPT"] = dual_hist["IPOPT"]
