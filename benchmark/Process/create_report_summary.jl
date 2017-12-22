include("../benchmark.jl")
include("create_report.jl")
overlapping_results = get_CUTEst_results()

#############
############# SUMMARY TABLE
#############

problem_list = collect(keys(first(overlapping_results)[2]))
method_list = collect(keys(overlapping_results))

different_status_problems = Array{String,1}()

begin
  for problem_name in problem_list

      previous_status = :none
      same_status = false

      for method_name in method_list
        info = overlapping_results[method_name][problem_name]

        if info.status == previous_status
          same_status = true
        end
        previous_status = info.status
      end

      if !same_status
        push!(different_status_problems,problem_name)
      end
  end
end

shared_fails = shared_failures(overlapping_results)
fail_only = Dict()
for (method_name, sum_data) in overlapping_results
    fail_only[method_name] = tot_failures(sum_data) - shared_fails
    println("total failures of $method_name is ", tot_failures(sum_data))
end

different_status_results = restrict_problems(overlapping_results,different_status_problems)


## create data frame
using DataFrames
df = DataFrame(
column_name = [
"failure",
"infeasible",
"unbounded",
#"worse obj",
#"better obj",
"KKT"
],
ipopt = [
fail_only["ipopt"],
tot(different_status_results["ipopt"],[:primal_infeasible]),
tot(different_status_results["ipopt"],[:dual_infeasible]),
#winners["one phase"],
#winners["ipopt"],
tot(different_status_results["ipopt"],[:optimal])
],
one_phase = [
fail_only["one phase"],
tot(different_status_results["one phase"],[:primal_infeasible]),
tot(different_status_results["one phase"],[:dual_infeasible]),
#winners["ipopt"],
#winners["one phase"],
tot(different_status_results["one phase"],[:optimal])
]
)


using CSV
CSV.write("$folder/summary.csv",df)


## bar plot of data frame
#Pkg.clone("https://github.com/JuliaPlots/StatPlots.jl.git")
using StatPlots
red2 = RGBA(1.0,0.0,0.0,0.6) #[red2,red2]
green2 = RGBA(1.0,65.0/256.0,256.0/256.0,63.0/256.0)
colour = [[:red,:red] [:orange,:orange] [:blue,:blue] [:green,:green]]
label = ["failure" "infeasible" "unbounded" "KKT"] #"worse obj" "better obj"
arr = Array(df)
data = convert(Array{Float64,2},arr[:,2:3])'
space = " "^50

for i = 1:length(label)
    groupedbar(data[:,1:i], bar_position = :stack,
    bar_width=0.9,
    colour=colour[:,1:i],
    label=label[:,1:i],
    xticks=[],
    xlabel="ipopt $space one phase",
    ylabel="number of problems",
    ylims=(0,80)
    )
    savefig("$folder/bar_$i.pdf")
end
