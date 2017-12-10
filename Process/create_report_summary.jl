include("../include.jl")
include("create_report.jl")
overlapping_results = get_CUTEst_results()

f_TOL = 1e-2

problem_list = collect(keys(first(overlapping_results)[2]))
method_list = collect(keys(overlapping_results))

winners = Dict()
winners["tie"] = 0
winners["fail"] = 0

different_status_problems = Array{String,1}()

begin
  for method_name in method_list
      print(pd(""), pd(method_name))
      winners[method_name] = 0
  end
  print("\n")

  for problem_name in problem_list

      previous_status = :none
      same_status = false

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

        if info.status == previous_status
          same_status = true
        end
        previous_status = info.status

        if info.status == :optimal
          if con_vio > 1e-6
            both_kkt = false
            println("CON VIOLATION!!")
          else
            @assert(con_vio < 1e-6)
            ftol_scaled = f_TOL * (1.0 + abs(fbest))
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

      if !same_status
        push!(different_status_problems,problem_name)
      end
  end
end

shared_fails = shared_failures(overlapping_results)
fail_only = Dict()
for (method_name, sum_data) in overlapping_results
    fail_only[method_name] = tot_failures(sum_data) - shared_fails
end

different_status_results = restrict_problems(overlapping_results,different_status_problems)


## create data frame
using DataFrames
df = DataFrame(
column_name = [
"failure",
"infeasible",
"unbounded",
"worse obj",
"better obj",
"KKT"
],
ipopt = [
fail_only["ipopt"],
tot(different_status_results["ipopt"],[:primal_infeasible]),
tot(different_status_results["ipopt"],[:dual_infeasible]),
winners["one phase"],
winners["ipopt"],
tot(different_status_results["ipopt"],[:optimal])
],
one_phase = [
fail_only["one phase"],
tot(different_status_results["one phase"],[:primal_infeasible]),
tot(different_status_results["one phase"],[:dual_infeasible]),
winners["ipopt"],
winners["one phase"],
tot(different_status_results["one phase"],[:optimal])
]
)




## bar plot of data frame
#Pkg.clone("https://github.com/JuliaPlots/StatPlots.jl.git")
using StatPlots
red2 = RGBA(1.0,0.0,0.0,0.6) #[red2,red2]
green2 = RGBA(1.0,65.0/256.0,256.0/256.0,63.0/256.0)
colour = [[:red,:red] [:orange,:orange] [:purple,:purple] [:blue,:blue] [:green,:green] [green2, green2]]
label = ["failure" "infeasible" "unbounded" "worse obj" "better obj" "KKT"]
arr = Array(df)
data = convert(Array{Float64,2},arr[:,2:3])'
space = " "^50

for i = 1:6
    groupedbar(data[:,1:i], bar_position = :stack,
    bar_width=0.9,
    colour=colour[:,1:i],
    label=label[:,1:i],
    xticks=[],
    xlabel="ipopt $space one phase",
    ylabel="number of problems",
    ylims=(0,80)
    )
    savefig("output/bar_$i.pdf")
end
