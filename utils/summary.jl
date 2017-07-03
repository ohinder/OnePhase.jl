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


function write_pars(stream::IOStream, par::Class_parameters)
    write(stream, "PAR \t\t\t\t\t\t VALUE \n")

    for fieldname in fieldnames(par)
        write(stream, "$(pd(fieldname,40)) \t $(getfield(par, fieldname)) \n")
    end
end

function write_summary(stream::IOStream, summary::Dict{String, problem_summary})
    write(stream, "problem_name \t status \t it count \t total_time \n");
    for (problem_name, info) in summary
      write(stream, "$(pd(problem_name)) \t $(pd(info.status)) \t $(pd(info.it_count)) \t $(rd(info.total_time))\n");
    end
end

function set_info_me!(info::problem_summary, hist)
    info.it_count = hist[end].t
    info.fval = hist[end].fval;
    info.con_vio = hist[end].con_vio;
    info.dual_feas = hist[end].norm_grad_lag;
    info.comp = hist[end].comp;
end
