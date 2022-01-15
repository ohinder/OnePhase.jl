# code for taking two tables of results and doing comparisions
using CSV, DataFrames,  PyPlot

function df_append_columns!(df::DataFrame,append::String;ignore=Array{Symbol,1}())
    for old_name in names(df)
        new_name = String(old_name) * append
        if !(old_name in ignore)
            rename!(df,old_name => Symbol(new_name))
        end
    end
end

function status_table(df_combine)
    # create a CSV with the statuses of each solver
end

function ratio_to_fastest(mat::Matrix)
    # loop through all columns of the dataframe and pick
    #@assert(mat .> 0.0)
    ratios = zeros(size(mat))
    for i = 1:size(mat,1)
        fastest_val = 10.0^10.0
        for j = 1:size(mat,2)
            if isnan(fastest_val)
                warn("NaN")
            else
                #@show mat[i,j]
                fastest_val = min(fastest_val,mat[i,j])
            end
        end

        for j = 1:size(mat,2)
            if mat[i,j] > 0.0 && fastest_val > 0.0
                ratios[i,j] = mat[i,j] / fastest_val
            end
        end
    end

    return ratios
end

function sort_mat!(mat::Matrix)
    for j = 1:size(mat,2)
        mat[:,j] = sort(mat[:,j])
    end
end

function plot_ratios!(ratios::Matrix;max_ratio=128.0,labels=Array{String,1}(),colors=Array{String,1}(),lines=Array{String,1}())
    sort_mat!(ratios)
    n = size(ratios,1)

    for j = 1:size(ratios,2)
        vals = ratios[:,j]
        only_plot = vals[vals .> 1.0];
        theta = (length(vals) - length(only_plot)) / length(vals)
        y = linspace(theta,1,length(only_plot))
        if length(labels) > 0
            label = labels[j]
        else
            label = ""
        end
        semilogx(only_plot,y,basex=2,label=label)
    end

    ylim([0.0,1.0])
    xlim([1.0,max_ratio])
    if length(labels) > 0
        legend()
    end
end

function geomean(vals::Vector)
    n = length(vals)
    return exp((1/n) * sum(log.(vals)))
end

function summary(df::DataFrame)
    data = Dict{Symbol,Array{Union{String,Float64}}}()
    data[:stat] = ["median","mean","geomean"];

    for h in names(df)
        data[h]=[median(df[:,h]),mean(df[:,h]),geomean(df[:,h])]
    end

    df_new = DataFrame(data)

    return df_new
end

function count(vec::Vector)
    count_dic = Dict()
    for el in vec
        if !haskey(count_dic,el)
            count_dic[el] = 1
        else
            count_dic[el] += 1
        end
    end
    return count_dic
end

function success(vec::Vector)
    return (vec .== "Optimal") .| (vec .== "dual_infeasible") .| (vec .== "primal_infeasible") .| (vec .== "Optimal_SS")
end

function compare_status(vec1::Vector,vec2::Vector)

end

test_name = "CUTEst/test-July-1"
RESULTS_DIR = "../benchmark-2/results/$test_name"

df_Ipopt = readtable("$RESULTS_DIR/Ipopt_summary.csv");
df_append_columns!(df_Ipopt,"_Ipopt",ignore=[:name])
df_OnePhase = readtable("$RESULTS_DIR/OnePhase_summary.csv");
df_append_columns!(df_OnePhase,"_OnePhase",ignore=[:name])

df_combine = join(df_Ipopt, df_OnePhase, on = :name) #, makeunique=true)
CSV.write("$RESULTS_DIR/combined.csv",df_combine)
df_combine = df_combine[df_combine[:status_Ipopt] .!= "TOO_FEW_DEGREES_OF_FREEDOM",:]
df_combine[df_combine[:status_Ipopt] .== "Optimal_SS",:status_Ipopt] = :Optimal;

count(df_combine[:status_Ipopt])
count(df_combine[:status_OnePhase])
count(df_combine[:status_Ipopt] .* df_combine[:status_OnePhase])

both_optimal_bools = (df_combine[:status_Ipopt].=="Optimal") .& (df_combine[:status_OnePhase].=="Optimal")
df_both_optimal = df_combine[both_optimal_bools,:];

fac_ratios = ratio_to_fastest([df_both_optimal[:n_fac_Ipopt] df_both_optimal[:n_fac_OnePhase]]);
plot_ratios!(fac_ratios,labels=["#fac Ipopt","#fac One Phase"])

fval_ratios = ratio_to_fastest([df_both_optimal[:n_fval_Ipopt] df_both_optimal[:n_fval_OnePhase]]);
plot_ratios!(fval_ratios,labels=["#fval Ipopt","#fval One Phase"])

n_hess_ratios = ratio_to_fastest([df_both_optimal[:n_Hess_Ipopt] df_both_optimal[:n_Hess_OnePhase]]);
plot_ratios!(n_hess_ratios,labels=["#hess Ipopt","#hess One Phase"])


df_n = df_both_optimal[:,[:iter_Ipopt,:iter_OnePhase,:n_fac_Ipopt,:n_fac_OnePhase,:n_fval_Ipopt,:n_fval_OnePhase,:n_Hess_Ipopt,:n_Hess_OnePhase]];
summary(df_n)

##################

function test_ratio_to_fastest()
    mat = rand(10,3)
    ratios = ratio_to_fastest(mat)
end



ratio_to_fastest([df_combine[:iter_Ipopt] df_combine[:iter_OnePhase]])


df_combine[:iter] ./ df_combine[:iter_1]

df_opt = df_combine[df_combine[:status] .== "Optimal",:]

it_ratio = log.(df_opt[:iter]) - log.(df_opt[:iter_1])
fac_ratio = log.(df_opt[:n_fac]) - log.(df_opt[:n_fac_1])
f_ratio = log.(df_opt[:n_fval]) - log.(df_opt[:n_fval_1])
c_ratio = log.(df_opt[:n_cons]) - log.(df_opt[:n_cons_1])

median(exp(it_ratio))
median(exp(fac_ratio))
median(exp(f_ratio))
median(exp(c_ratio))



sum(df_combine[:status] .== "Optimal")
sum(df_combine[:status_1] .== "Optimal")
sum(df_combine[:status] .== "primal_infeasible")
sum(df_combine[:status_1] .== "primal_infeasible")
sum(df_combine[:status] .== "dual_infeasible")
sum(df_combine[:status_1] .== "dual_infeasible")
