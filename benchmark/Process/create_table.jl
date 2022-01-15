include("../benchmark.jl")
include("create_report.jl")

results = get_CUTEst_results()
#["one_phase"] = load("../results/one_phase/Dec10/summary.jld", "summary")


using DataFrames,CSV

df_results = Dict()
for (method_name,method_results) in results
  df_results[method_name] = DataFrame(name=Array{String,1}(),it=Array{Int64,1}(),time=Array{Float64,1}(),fval=Array{Float64,1}(),con=Array{Float64,1}(),status=[])
  for (name,info) in method_results
    push!(df_results[method_name],[name,info.it_count,info.total_time,info.fval,info.con_vio,string(info.status)])
  end
  CSV.write("$folder/table_$(folder)_$(method_name).csv", df_results[method_name])
end



###
### ITER/TIMING BREAK DOWN
###

df_summary = DataFrame(name=[],median_it=[],median_time=[],total_time=[])
for (method_name,method_results) in results
    push!(df_summary,[method_name,median(df_results[method_name][:it]),median(df_results[method_name][:time]),sum(df_results[method_name][:time])])
end

df_summary

CSV.write("$folder/runtime_summary.csv", df_summary)

using FreqTables

##
## FAILURE REASONS
##
freqtable(df_results["ipopt"][:status])
freqtable(df_results["one phase"][:status])
