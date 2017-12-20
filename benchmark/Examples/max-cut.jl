using JuMP
include("../benchmark.jl")

function max_cut_model(n::Int64,r::Int64,edge_list::Array{Tuple{Int64,Int64},1})
    #edges = [(i,i+1) for i = 1:(n-1)]
    m = Model()

    #start_mat = randn(n,r)
    #
    #for i = 1:n
    #  start_mat[i,:] = start_mat[i,:] / norm(start_mat[i,:],2)
    #end

    @variable(m, x[i=1:n,k=1:r],start=randn() / sqrt(r))
    nedges = length(edge_list)
    @NLobjective(m, Min, sum(x[i,k] * x[j,k] for k = 1:r for (i,j) in edge_list) / 2.0 - nedges / 2.0)
    for i = 1:n
      @NLconstraint(m, sum(x[i,k]^2 for k=1:r) == 1)
    end

    return m
    #
end
edge_list = Array{Tuple{Int64,Int64},1}()
n = 200
for sample = 1:(20 * n)
  edge = (rand(1:n),rand(1:n))
  if edge[1] < edge[2] && !(edge in edge_list)
    push!(edge_list,edge)
  end
end

m = max_cut_model(n,6,edge_list)

pars = Class_parameters()
#pars.aggressive_dual_threshold = 1e-4
#pars.primal_bounds_dual_feas = true
iter, status, hist, t, err = one_phase_solve(m,pars);

using Ipopt
nlp_raw = MathProgNLPModel(m);
IpoptSolve(nlp_raw)
