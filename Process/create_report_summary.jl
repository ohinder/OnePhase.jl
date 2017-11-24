include("../include.jl")
include("create_report.jl")
results = get_CUTEst_results()

f_TOL = 1e-1

problem_list = collect(keys(first(overlapping_opt_results)[2]))
method_list = collect(keys(overlapping_opt_results))

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
        info = overlapping_opt_results[method_name][problem_name]
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
  end
end

shared_fails = shared_failures(overlapping_results)
fail_only = Dict()
for (method_name, sum_data) in overlapping_results
    fail_only[method_name] = tot_failures(sum_data) - shared_fails
end

using DataFrames
df = DataFrame(
column_name = [
"shared failures",
"only failure",
"worse KKT",
"infeasible",
"unbounded"
],
ipopt = [
shared_fails,
fail_only["ipopt"],
winners["one phase"],
tot(overlapping_results["ipopt"],[:primal_infeasible]),
tot(overlapping_results["ipopt"],[:dual_infeasible]),
],
one_phase = [
shared_fails,
fail_only["one phase"],
winners["ipopt"],
tot(overlapping_results["one phase"],[:primal_infeasible]),
tot(overlapping_results["one phase"],[:dual_infeasible]),
]
)





#Pkg.clone("https://github.com/JuliaPlots/StatPlots.jl.git")

using StatPlots
red2 = RGBA(1.0,0.0,0.0,0.6)
colour = [[:red,:red] [red2,red2] [:orange,:orange] [:green,:green] [:blue,:blue]]
label = ["both fail" "one fails" "infeasible" "unbounded" "worse KKT"]
groupedbar(rand(2,5), bar_position = :stack, bar_width=0.7, colour=colour, label=label)
