type problem_summary
    status::Symbol
    it_count::Int64
    total_time::Float64
    fval::Float64
    con_vio::Float64
    dual_feas::Float64
    comp::Float64

    function problem_summary()
        return new(:Blank,-2,NaN,NaN,NaN,NaN)
    end
end

abstract abstract_alg_history

#= 31st May
type problem_summary
    status::Symbol
    it_count::Int64
    total_time::Float64

    function problem_summary()
        return new()
    end
end
=#

#::IOStream
function write_pars(stream, par::Class_parameters)
    write(stream, "PAR \t\t\t\t\t\t VALUE \n")

    for fieldname in fieldnames(par)
        fieldval = getfield(par, fieldname)
        if isa(fieldval,abstract_pars)
          write(stream, "$(pd(fieldname,40)) \n")
          for subfieldname in fieldnames(fieldval)
            subfieldval = getfield(fieldval, subfieldname)
            write(stream, "\t $(pd(subfieldname,40)) \t $(subfieldval) \n")
          end
        else
          write(stream, "$(pd(fieldname,40)) \t $(fieldval) \n")
        end
    end
end

function write_summary(stream::IOStream, summary::Dict{String, problem_summary})
    write(stream, "problem_name \t status \t it count \t total_time \n");
    for (problem_name, info) in summary
      write(stream, "$(pd(problem_name)) \t $(pd(info.status)) \t $(pd(info.it_count)) \t $(rd(info.total_time))\n");
    end
end

function set_info_me!(info::problem_summary, hist, status::Symbol)
    info.it_count = hist[end].t
    info.fval = hist[end].fval;
    info.con_vio = hist[end].con_vio;
    info.dual_feas = hist[end].norm_grad_lag;
    info.comp = hist[end].comp;
    info.status = status
end

function set_cutest_info_ipopt!(info::problem_summary, ipopt_solver, nlp_raw::AbstractNLPModel, x::Array{Float64})
  num_vars = length(nlp_raw.meta.lvar)
  x_true = x[1:num_vars]

  info.fval = obj(nlp_raw, x_true)
  a = cons(nlp_raw, x_true);
  info.con_vio = max(0.0, maximum(nlp_raw.meta.lvar - x_true), maximum(x_true - nlp_raw.meta.uvar), maximum(nlp_raw.meta.lcon - a), maximum(a - nlp_raw.meta.ucon))

  info.dual_feas = norm(grad(nlp_raw, x_true) + jac(nlp_raw, x_true)' * ipopt_solver.mult_g + ipopt_solver.mult_x_U - ipopt_solver.mult_x_L, Inf);
  info.comp = maximum(ipopt_solver.mult_x_U .* min(1e16, abs(nlp_raw.meta.uvar - x_true) ) ) + maximum( ipopt_solver.mult_x_L .* min(1e16, abs(x_true - nlp_raw.meta.lvar)) )
end
