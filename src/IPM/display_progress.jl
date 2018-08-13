@compat abstract type abstract_alg_history end

export abstract_alg_history, get_col, generic_alg_history, alg_history2

using NLPModels, JuMP

struct alg_history2 <: abstract_alg_history
    t::Int64
    step_type::String
    ls_info_num_steps::Int64
    ls_info_step_size_P::Float64
    ls_info_step_size_D::Float64
    dir_x_norm::Float64
    dir_y_norm::Float64
    dir_s_norm::Float64
    kkt_ratio::Float64
    mu::Float64
    fval::Float64
    dual_scaled::Float64
    norm_grad_lag::Float64
    primal_residual::Float64
    con_vio::Float64
    comp::Float64
    comp_ratio::Float64
    dot_sy::Float64
    farkas::Float64
    delta::Float64
    eval_merit_function::Float64
    eval_phi::Float64
    eval_grad_phi::Float64
    y_norm::Float64
    x_norm::Float64
    tot_num_fac::Int64
    num_fac_inertia::Int64
    strict_comp::Float64

    #=function alg_history()
      return new()
  end=#
end

function major_its_only(hist::Array{alg_history2,1})
    # reduce the history list so it only includes major iterations
    its = get_col(hist,:t);
    keep_these = Array{Int64,1}()
    for i = 1:(length(its)-1)
        if its[i+1] != its[i]
            push!(keep_these,i)
        end
    end
    push!(keep_these,length(its))
    return hist[keep_these]
end

type generic_alg_history <: abstract_alg_history
  t::Int64
  #mu::Float64
  fval::Float64
  #dual_scaled::Float64
  norm_grad_lag::Float64
  comp::Float64
  #comp_ratio::Float64
  #dot_sy::Float64
  #primal_residual::Float64
  #farkas::Float64
  #delta::Float64
  #eval_merit_function::Float64
  #eval_phi::Float64
  #eval_grad_phi::Float64
  primal_residual::Float64
  con_vio::Float64
  y_norm::Float64
  x_norm::Float64
  strict_comp::Float64
end

function get_col(arr::Array,colname::Symbol)
    col = nothing
    for i = 1:length(arr)
        val = getfield(arr[i],colname)
        if i == 1
            col = Array{typeof(val),1}();
        end
        push!(col,val)
    end
    return col
end

function head_progress()
  println( pd("it",3), pd("step", 6), pd("eta",8), pd("α_P"), pd("α_D"), pd("ls",3), pd("|dx|"), pd("|dy|"), pd("N err"), "|", pd("mu"), pd("dual"), pd("primal"), pd("cmp scaled"), pd("infeas?"), "|", pd("delta"), pd("#fac",6), pd("|x|"), pd("|y|"), pd("∇phi"), pd("phi") )
end

function record_progress_first_it!(progress::Array{alg_history2,1}, iter::Class_iterate, kkt_solver::abstract_KKT_system_solver, par::Class_parameters, display::Bool)
    record_progress!(progress, 0, "it0", iter, kkt_solver, Blank_ls_info(), Class_reduction_factors(),  0, 0, par, display)
end


function record_progress!(progress::Array{alg_history2,1}, t::Int64, step_type::String, iter::Class_iterate, kss::abstract_KKT_system_solver,  ls_info::abstract_ls_info, eta::Class_reduction_factors, num_fac_inertia::Int64, tot_num_fac::Int64, par::Class_parameters, display::Bool)
    dir_y_norm = norm(kss.dir.y,Inf);
    dir_x_norm = norm(kss.dir.x,Inf);
    dir_s_norm = norm(kss.dir.s,Inf);
    ls_info_num_steps = ls_info.num_steps
    ls_info_step_size_P = ls_info.step_size_P
    ls_info_step_size_D = ls_info.step_size_D
    kkt_err = kss.kkt_err_norm;
    mu = get_mu(iter)
    fval = get_fval(iter)
    dot_sy = dot(iter.point.s,iter.point.y)/length(iter.point.s)
    norm_grad_lag = norm(eval_grad_lag(iter, 0.0),Inf)
    norm_grad_lag_mod = norm(eval_grad_lag(iter, iter.point.mu),Inf)
    dual_scaled = scaled_dual_feas(iter, par)
    primal_residual = norm(get_primal_res(iter),Inf)
    con_vio = get_max_vio(iter)
    val_comp = norm(comp(iter),Inf)
    comp_ratio = comp_ratio_max(iter) * par.ls.comp_feas_agg #norm(comp(iter),Inf)
    farkas = eval_farkas(iter)
    delta = get_delta(iter)
    val_merit_function = eval_merit_function(iter, par)
    val_phi = eval_phi(iter,iter.point.mu)
    val_grad_phi = norm(eval_grad_phi(iter, iter.point.mu),Inf)
    y_norm = norm(get_y(iter),Inf)
    x_norm = norm(get_x(iter),Inf)
    tot_num_fac = tot_num_fac
    num_fac_inertia = num_fac_inertia
    kkt_ratio = kss.kkt_err_norm.ratio
    strict_comp = minimum(iter.point.s+iter.point.y)


    hist = alg_history2(t::Int64,
    step_type::String,
    #eta::Class_reduction_factors,
    ls_info_num_steps::Int64,
    ls_info_step_size_P::Float64,
    ls_info_step_size_D::Float64,
    dir_x_norm::Float64,
    dir_y_norm::Float64,
    dir_s_norm::Float64,
    kkt_ratio::Float64,
    mu::Float64,
    fval::Float64,
    dual_scaled::Float64,
    norm_grad_lag::Float64,
    primal_residual::Float64,
    con_vio::Float64,
    val_comp::Float64,
    comp_ratio::Float64,
    dot_sy::Float64,
    farkas::Float64,
    delta::Float64,
    val_merit_function::Float64,
    val_phi::Float64,
    val_grad_phi::Float64,
    y_norm::Float64,
    x_norm::Float64,
    tot_num_fac::Int64,
    num_fac_inertia::Int64,
    strict_comp::Float64)

    push!(progress, hist)

    if display
      println(pd(t,5), pd(step_type[1], 3), rd(mu),  rd(ls_info_step_size_P),  rd(ls_info_step_size_D), pd(ls_info_num_steps,2), rd(dir_x_norm), rd(dir_y_norm), rd(kkt_ratio), "|", rd(mu), rd(dual_scaled), rd(primal_residual), rd(comp_ratio), rd(farkas), "|", rd(delta), pd(num_fac_inertia,3), pd(tot_num_fac,3), rd(x_norm), rd(y_norm), rd(val_grad_phi), rd(val_merit_function, 5))
    end
end
