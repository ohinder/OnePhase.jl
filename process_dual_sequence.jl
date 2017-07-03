include("include.jl")

summary_netlib = Dict{String, Dict{String,problem_summary}}()
summary_netlib["netlib/well-behaved"] = load("results/netlib/well-behaved/summary.jld", "summary")
#summary_netlib["netlib/fast-primal"] = load("results/netlib/fast-primal/summary.jld", "summary")
summary_netlib["netlib/ipopt-no-perturb"] = convert_JuMP(load("results/netlib/ipopt-no-perturb/summary.jld", "summary"))
summary_netlib["netlib/ipopt-perturb"] = convert_JuMP(load("results/netlib/ipopt-perturb/summary.jld", "summary"))


for (method_name, sum_data) in summary_netlib
    println(method_name, "=" ,tot_failures(sum_data), " out of ", length(sum_data))
end


problem_list, method_list = get_lists(summary_netlib)
all_opt_problems = get_all_optimal(summary_netlib, problem_list)

function compute_dual_maxes(summary, folder_name)
    dual_maxes = Dict{String,Dict{String,Float64}}() #Array{Float64,1}
    for (method_name, sum_data) in summary
      dual_maxes[method_name] = Dict{String,Float64}()
      for (problem_name, info) in sum_data
          file = "$folder_name/$method_name/jld/$problem_name.jld"
          hist = load(file, "history")
          dual_hist = compute_dual_hist(hist)
          dual_maxes[method_name][problem_name] = maximum(dual_hist)
      end
    end

    return dual_maxes
end

linestyle_dict = Dict(
  "netlib/well-behaved" => "--",
  "netlib/ipopt-no-perturb" => "-",
  "netlib/fast-primal" => "-.",
  "netlib/ipopt-perturb" => ":"
)

label_dic = Dict(
  "netlib/well-behaved" => "Well behaved IPM",
  "netlib/ipopt-no-perturb" => "IPOPT w/o perturb",
  "netlib/fast-primal" => "IPM reduces primal fast",
  "netlib/ipopt-perturb" => "IPOPT w. perturb"
)

dual_maxes_all = compute_dual_maxes(summary_netlib, "results")

dual_maxes_restricted = Dict{String,Array{Float64,1}}()
for method_name in method_list
  dual_maxes_restricted[method_name] = Array{Float64,1}()
  for problem_name in all_opt_problems
   push!(dual_maxes_restricted[method_name], dual_maxes_all[method_name][problem_name])
  end
end


using PyPlot

min_y = 1.0

y_vals = collect(1:length(all_opt_problems)) / length(all_opt_problems)
for (method_name, val) in dual_maxes_restricted
  semilogy(y_vals, sort(dual_maxes_restricted[method_name]),linestyle=linestyle_dict[method_name], label=label_dic[method_name], color="black", basey=10)
end

xlabel("proportion of problems")
ylabel("maximum dual value over all iterates")

legend()


#df = DataFrame(problem_name = problem_list, feasible11 = status11, problem_name2 = names2, feasible21 = status21 )


# IPOPT 11/29 failures
# us 4/29 failures




con_vio = Dict{String,Array{Float64,1}}()
for folder_name in method_list
  con_vio[folder_name] = Array{Float64,1}()
  for problem_name in problem_list
      results = load("results/$folder_name/jld/$problem_name.jld", "history")
      push!(con_vio[folder_name], results[end].con_vio)
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

####################################
### PLOT ONE PARTICULAR PROBLEM
####################################

field_list = Dict("max dual variable" => :y_norm, "primal feasibility" => :con_vio, "dual feasibility" => :norm_grad_lag, "complementarity" => :comp);
line_style_list = Dict("max dual variable" => "-", "primal feasibility" => "--", "dual feasibility" => "-.", "complementarity" => ":");

method_name = "well-behaved"
problem_name = "ADLITTLE"
file_wb = "results/netlib/$method_name/jld/$problem_name.jld"
hist_wb = load(file_wb, "history")
hist_by_it_wb = summarize_by_iteration(hist_wb)
data_wb = aggregate_summary_by_it_for_plots(hist_by_it_wb, field_list)


method_name = "ipopt-no-perturb"
file_ipnp = "results/netlib/$method_name/jld/$problem_name.jld"
hist_ipnp = load(file_ipnp, "history")
hist_by_it_ipnp = summarize_by_iteration(hist_ipnp)
data_ipnp = aggregate_summary_by_it_for_plots(hist_by_it_ipnp, field_list)

method_name = "ipopt-perturb"
file_ipp = "results/netlib/$method_name/jld/$problem_name.jld"
hist_ipp = load(file_ipp, "history")
hist_by_it_ipp = summarize_by_iteration(hist_ipp)
data_ipp = aggregate_summary_by_it_for_plots(hist_by_it_ipp, field_list)


##
## DO THE PLOTTING
##

function plot_single_example_compare_iterates(data, field_list, line_style_list)
  for (label, field) in field_list
    semilogy(data[field],label=label,linestyle=line_style_list[label], color="black")
  end
end

ylims = [1e-9, 1e12]

#plt[:figlegend]( ["-"], ["test"], loc = (0.5, 0), ncol=5 )
#plt[:figlegend]( ["-"], ["test"], loc = (0.0, 0), ncol=1 )




subplot(131)
plot_single_example_compare_iterates(data_wb, field_list, line_style_list)

title("A well-behaved IPM")
ylabel("L-infinity norm")
xlabel("iterations")
ax = gca()
ax[:set_ylim](ylims)

subplot(132)
plot_single_example_compare_iterates(data_ipnp, field_list, line_style_list)

title("IPOPT w/o perturbation")
#ylabel("L-infinity norm")
xlabel("iterations")
ax = gca()
ax[:set_ylim](ylims)
ax[:set_yticklabels]([])

subplot(133)
plot_single_example_compare_iterates(data_ipp, field_list, line_style_list)

title("IPOPT w. perturbation")
xlabel("iterations")
ax = gca()
ax[:set_ylim](ylims)
ax[:set_yticklabels]([])

legend()
ax[:legend]()[:draggable]()

#suptitle("Comparsion of algorithms on $problem_name")
#ttl = ax[:title]
#ttl[:set_position]([.5, 1.05])
#plt[:tight_layout]()

#ax[:set_position]([.5, 1.05])
