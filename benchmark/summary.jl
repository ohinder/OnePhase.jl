type problem_summary # for backward compatibility
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

function convert_problem_summary(ps::problem_summary)
    ps_new = problem_summary2()
    ps_new.status = ps.status
    ps_new.it_count = ps.it_count
    ps_new.total_time = ps.total_time
    ps_new.fval = ps.fval
    ps_new.con_vio = ps.con_vio
    ps_new.dual_feas = ps.dual_feas
    ps_new.comp = ps.comp
    ps_new.dual_max = NaN

    return ps_new
end

function convert_problem_summary(ls::Dict{String,problem_summary})
    new_ls = Dict{String,problem_summary2}()
    for (problem_name,info) in ls
      new_ls[problem_name] = info
    end

    return new_ls
end

type problem_summary2
    status::Symbol
    it_count::Int64
    total_time::Float64
    fval::Float64
    con_vio::Float64
    dual_feas::Float64
    comp::Float64
    dual_max::Float64

    function problem_summary2()
        return new(:Blank,-2,NaN,NaN,NaN,NaN,NaN)
    end
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

function write_summary(stream::IOStream, summary::Dict{String, problem_summary2})
    write(stream, "problem_name \t status \t it count \t total_time \n");
    for (problem_name, info) in summary
      write(stream, "$(pd(problem_name)) \t $(pd(info.status)) \t $(pd(info.it_count)) \t $(rd(info.total_time))\n");
    end
end

function write_summary(file_name::String, summary::Dict{String, problem_summary2})
    summary_file = open(file_name, "w")
    write_summary(summary_file, summary)
    close(summary_file)
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

function set_cutest_info_ipopt!(info::problem_summary2, ipopt_solver, nlp_raw::AbstractNLPModel, x::Array{Float64})
  num_vars = length(nlp_raw.meta.lvar)
  x_true = x[1:num_vars]

  info.fval = obj(nlp_raw, x_true)
  a = cons(nlp_raw, x_true);
  info.con_vio = max(0.0, maximum(nlp_raw.meta.lvar - x_true), maximum(x_true - nlp_raw.meta.uvar), maximum(nlp_raw.meta.lcon - a), maximum(a - nlp_raw.meta.ucon))

  info.dual_feas = norm(grad(nlp_raw, x_true) + jac(nlp_raw, x_true)' * ipopt_solver.mult_g + ipopt_solver.mult_x_U - ipopt_solver.mult_x_L, Inf);
  info.comp = maximum(ipopt_solver.mult_x_U .* min(1e16, abs(nlp_raw.meta.uvar - x_true) ) ) + maximum( ipopt_solver.mult_x_L .* min(1e16, abs(x_true - nlp_raw.meta.lvar)) )
  info.dual_max = max(maximum(ipopt_solver.mult_x_U),maximum(ipopt_solver.mult_x_L), abs(ipopt_solver.mult_g))
end
