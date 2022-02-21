include("include.jl")
include("read_lp.jl")

#println("GUROBI NOT WORKING")
#using Gurobi

### SHOULD REALLY BE LOOKING AT DUAL
###


#for tol in [1e-8]

function lp_feasible(tol, problem_list)
  status_list = Array{String,1}()
  for problem_name in problem_list
      mod = read_lp_into_JuMP(problem_name, -tol, false)
      setsolver(mod, GurobiSolver(FeasibilityTol=1e-9,BarHomogeneous=1))
      JuMP.build(mod)
      status =solve(mod)
      string_status = (status == :Optimal) ? "true" : "false"
      push!(status_list, string_status)
  end

  return status_list
end
tol_set = [10.0, 1e-4, 1e-6, 1e-8]
status_list = Dict{Float64,Array{String,1}}()
problem_list = load_netlib(Int(1e5))

for tol in tol_set
  status_list[tol] = lp_feasible(tol, problem_list)
end

#end
#using CSV
#CSV.write("seq-ipm-figure/table.csv",results_table)
# Pkg.clone("https://github.com/farrellm/YTables.jl.git")
using DataFrames, YTables

@show interior, no_interior

status_list = results[tol]
divide = Int(length(status_list)/2)

#tol1 =

tol1 = 1e-4
tol2 = 1e-6
tol3 = 1e-8
status11 = status_list[tol1][1:divide]
status12 = status_list[tol2][1:divide]
names1 = problem_list[1:divide]
status21 = status_list[tol1][(divide+1):end]
status22 = status_list[tol2][(divide+1):end]
names2 = problem_list[(divide+1):end]

df = DataFrame(problem_name1 = names1, feasible11 = status11, problem_name2 = names2, feasible21 = status21 )

#writetable("seq-ipm-figures/table.csv", df)

latex(df)

interior = Dict()
no_interior = Dict()
for tol in tol_set
  interior[tol] = sum("true" .== status_list[tol])
  no_interior[tol] = sum("false" .== status_list[tol])
end
