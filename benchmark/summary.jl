using LinearAlgebra, SolverCore
mutable struct problem_summary # for backward compatibility
    status::Symbol
    it_count::Int64
    total_time::Float64
    fval::Float64
    con_vio::Float64
    dual_feas::Float64
    comp::Float64

    function problem_summary()
        return new(:Blank,-2,NaN,NaN,NaN)
    end
end
mutable struct problem_summary2
    status::Symbol
    it_count::Int64
    total_time::Float64
    fval::Float64
    con_vio::Float64
    dual_feas::Float64
    comp::Float64
    dual_max::Float64

    #size
    number_variables::Int64
    number_constraints::Int64
    #toal evaluations
    total_fval_evaluation::Int64
    total_grad_evaluation::Int64
    total_jac_evaluation::Int64
    total_cons_evaluation::Int64
    total_hess_evaluation::Int64

    function problem_summary2()
        return new(:Blank,-2,NaN,NaN,NaN,NaN,NaN,0,0,0,0,0,0,0)
    end
end

function cps(ps::problem_summary)
    ps_new = problem_summary2()
    ps_new.status = ps.status
    ps_new.it_count = ps.it_count
    ps_new.total_time = ps.total_time
    ps_new.fval = ps.fval
    ps_new.con_vio = ps.con_vio
    ps_new.dual_feas = ps.dual_feas
    ps_new.comp = ps.comp
    ps_new.dual_max = NaN

    ps_new.number_variables = ps.number_variables
    ps_new.number_variables = ps.number_variables
    #toal evaluations
    ps_new.total_fval_evaluation = ps.total_fval_evaluation
    ps_new.total_grad_evaluation = ps.total_grad_evaluation
    ps_new.total_jac_evaluation = ps.total_jac_evaluation
    ps_new.total_cons_evaluation = ps.total_cons_evaluation
    ps_new.total_hess_evaluation = ps.total_hess_evaluation

    return ps_new
end

function cps(ls::Dict{String,problem_summary})
    new_ls = Dict{String,problem_summary2}()
    for (problem_name,info) in ls
      new_ls[problem_name] = cps(info)
    end

    return new_ls
end

function cps(ls::Dict{String,problem_summary2})
    return ls
end


#= 31st May
type problem_summary2
    status::Symbol
    it_count::Int64
    total_time::Float64

    function problem_summary2()
        return new()
    end
end
=#

function write_summary(stream::IOStream, summary::Dict{String, problem_summary2})
    write(stream, "problem_name \t status \t it count \t total_time \n");
    for (problem_name, info) in summary
      # write(stream, "$(OnePhase.pd(problem_name)) \t $(OnePhase.pd(info.status)) \t $(OnePhase.pd(info.it_count)) \t $(OnePhase.rd(info.total_time))\n");
       write(stream, "$(OnePhase.pd(problem_name)) \t $(OnePhase.pd(info.status)) \t $(OnePhase.pd(info.it_count)) \t $(OnePhase.rd(info.total_time))\n");
    end
end

function write_summary(file_name::String, summary::Dict{String, problem_summary2})
    summary_file = open(file_name, "w")
    write_summary(summary_file, summary)
    close(summary_file)
end

function set_info_me!(info::problem_summary2, hist, status::Symbol, iter)
    set_info_me!(info, hist, status)
    info.total_fval_evaluation = iter.nlp.nlp.counters.neval_obj
    info.total_grad_evaluation = iter.nlp.nlp.counters.neval_grad
    info.total_jac_evaluation = iter.nlp.nlp.counters.neval_jac
    info.total_cons_evaluation =  iter.nlp.nlp.counters.neval_cons
    info.total_hess_evaluation = iter.nlp.nlp.counters.neval_hess
end

function set_info_me!(info::problem_summary2, hist, status::Symbol)
    info.it_count = hist[end].t
    info.fval = hist[end].fval;
    info.con_vio = hist[end].con_vio;
    info.dual_feas = hist[end].norm_grad_lag;
    info.comp = hist[end].comp;
    info.dual_max = hist[end].y_norm
    info.status = status
end

function set_cutest_info_ipopt!(info::problem_summary2, stats::SolverCore.GenericExecutionStats{Float64, Vector{Float64}}, nlp_raw::AbstractNLPModel, x::Array{Float64,1})
  num_vars = length(nlp_raw.meta.lvar)
  x_true = x[1:num_vars]

  info.fval = stats.objective
  info.total_fval_evaluation = stats.counters.neval_obj
  info.total_grad_evaluation = stats.counters.neval_grad
  info.total_jac_evaluation = stats.counters.neval_jac
  info.total_cons_evaluation =  stats.counters.neval_cons
  info.total_hess_evaluation = stats.counters.neval_hess

  a = cons(nlp_raw, x_true);
  info.con_vio = max(0.0, maximum(nlp_raw.meta.lvar - x_true), maximum(x_true - nlp_raw.meta.uvar), maximum(nlp_raw.meta.lcon - a), maximum(a - nlp_raw.meta.ucon))

  info.dual_feas = norm(grad(nlp_raw, x_true) + jac(nlp_raw, x_true)' * stats.multipliers + stats.multipliers_U - stats.multipliers_L, Inf);
  info.comp = maximum(stats.multipliers_U .* min.(1e16, abs.(nlp_raw.meta.uvar - x_true) ) ) + maximum( stats.multipliers_L .* min.(1e16, abs.(x_true - nlp_raw.meta.lvar)) )
  info.dual_max = max(maximum(stats.multipliers_U),maximum(stats.multipliers_L), maximum(abs.(stats.multipliers)))
end
