#=type Class_ls_info
    step_size_P::Float64
    step_size_D::Float64
    num_steps::Int64
    frac_progress::Float64

    function Class_ls_info()
        return new(NaN,NaN,0,0.0)
    end

    function Class_ls_info(step_size_P::Float64, step_size_D::Float64, num_steps::Int64, frac_progress::Float64)
        return new(step_size_P, step_size_D, num_steps, frac_progress)
    end
end=#

function scale_direction(dir::Class_point, step_size::Float64)
    new_dir = deepcopy(dir)
    shrink_direction!(new_dir, step_size)
    return new_dir
end

#=function max_step_LP(iter::Class_iterate, dir::Class_point, pars::Class_parameters)
    for i = 1:ncon(iter)

    end
end=#
function max_step_primal_dual(iter::Class_iterate, dir::Class_point, threshold::Float64)
    return simple_max_step([iter.point.s; iter.point.y], [dir.s; dir.y], threshold)
end

function max_step_primal(iter::Class_iterate, dir::Class_point, threshold::Float64)
    return simple_max_step(iter.point.s, dir.s, threshold)
end

function simple_max_step(val::Array{Float64,1}, dir::Array{Float64,1}, threshold::Float64)
    q = 1.0 / (1.0 - threshold)
    ratio = maximum( [1.0; -q * dir ./ val ] ) #; -q * dir.y ./ iter.point.y] );
    return 1.0 / ratio
end

abstract abstract_ls_info;

type Class_stable_ls <: abstract_ls_info
    step_size_P::Float64
    step_size_D::Float64
    num_steps::Int64
    frac_progress::Float64
    predict_red::Float64
    actual_red::Float64
    cur_merit::Float64

    function Class_stable_ls(iter::Class_iterate, dir::Class_point, pars::Class_parameters)
        this = new()
        this.predict_red = merit_function_predicted_reduction(iter, dir);
        this.cur_merit = eval_merit_function(iter, pars)

        return this
    end

    function Class_stable_ls(step_size_P::Float64,step_size_D::Float64,num_steps::Int64,frac_progress::Float64, predict_red::Float64)
        this = new(step_size_P, step_size_D, num_steps, frac_progress, predict_red)
        return this
    end
end

function Blank_ls_info()
    this = Class_stable_ls(0.0,0.0, 0, 0.0, 0.0)
    return this
end

type Class_agg_ls <: abstract_ls_info
    step_size_P::Float64
    step_size_D::Float64
    num_steps::Int64
    frac_progress::Float64
    predict_red::Float64
    actual_red::Float64

    function Class_agg_ls(iter::Class_iterate, dir::Class_point, pars::Class_parameters)
        this = new()
        this.predict_red = -1.0
        return this
    end
end

type Class_filter_ls <: abstract_ls_info
    step_size_P::Float64
    step_size_D::Float64
    num_steps::Int64
    frac_progress::Float64
    predict_red::Float64
    actual_red::Float64
    cur_merit::Float64

    function Class_filter_ls(iter::Class_iterate, dir::Class_point, pars::Class_parameters)
        this = new()
        this.predict_red = merit_function_predicted_reduction(iter, dir);
        this.cur_merit = eval_merit_function(iter, pars)

        return this
    end

    function Class_filter_ls(step_size_P::Float64,step_size_D::Float64,num_steps::Int64,frac_progress::Float64, predict_red::Float64)
        this = new(step_size_P, step_size_D, num_steps, frac_progress, predict_red)
        return this
    end
end

type Class_kkt_ls <: abstract_ls_info
    step_size_P::Float64
    step_size_D::Float64
    num_steps::Int64
    frac_progress::Float64
    predict_red::Float64
    actual_red::Float64
    cur_merit::Float64

    function Class_kkt_ls(iter::Class_iterate, dir::Class_point, pars::Class_parameters)
        this = new()
        this.predict_red = merit_function_predicted_reduction(iter, dir);
        this.cur_merit = eval_merit_function(iter, pars)

        return this
    end

    function Class_kkt_ls(step_size_P::Float64,step_size_D::Float64,num_steps::Int64,frac_progress::Float64, predict_red::Float64)
        this = new(step_size_P, step_size_D, num_steps, frac_progress, predict_red)
        return this
    end
end

type Class_filter
    fval::Float64
    scaled_kkt_err::Float64
    mu::Float64

    function Class_filter(iter::Class_iterate, pars::Class_parameters)
      this = new()
      this.fval = eval_phi(iter)
      this.scaled_kkt_err = eval_kkt_err(iter, pars)
      this.mu = iter.point.mu

      return this
    end
end

function add!(ar::Array{Class_filter,1}, iter::Class_iterate, pars::Class_parameters)
    push!( ar, Class_filter(iter, pars) )
end

function satisfies_filter!(ar::Array{Class_filter,1}, can::Class_iterate, step_size::Float64, pars::Class_parameters)
    p = Class_filter(can, pars)

    for i = 1:length(ar)
        #tol = 1e-3 * max(1.0, p.scaled_kkt_err)
        expected_reduction = pars.kkt_reduction_factor
        #p.fval < ar[i].fval - (p.scaled_kkt_err)^2 ||
        if !( p.mu < ar[i].mu || p.scaled_kkt_err / ar[i].scaled_kkt_err < (1.0 - pars.kkt_reduction_factor * step_size))
            return false
        end
    end

    return true
end

function accept_func_stable!(accept::abstract_ls_info, iter::Class_iterate, candidate::Class_iterate, pars::Class_parameters)
  if accept.predict_red < 0.0
    old_merit = eval_merit_function(iter, pars)
    new_merit = eval_merit_function(candidate, pars)
    accept.actual_red = new_merit - old_merit
    accept.frac_progress = accept.actual_red / accept.predict_red
    if accept.frac_progress > pars.predict_reduction_factor * accept.step_size_P
        #@show accept.frac_progress, accept.actual_red
        return :success
    else
        return :not_enough_progress
    end
  else
      return :not_enough_progress
  end
end

function accept_func!(accept::Class_stable_ls, iter::Class_iterate, candidate::Class_iterate, filter::Array{Class_filter,1}, pars::Class_parameters)
    status = accept_func_stable!(accept, iter, candidate, pars)
    return status
end

function accept_func_kkt!(accept::abstract_ls_info, iter::Class_iterate, candidate::Class_iterate, filter::Array{Class_filter,1}, pars::Class_parameters)

    if satisfies_filter!(filter, candidate, accept.step_size_P, pars)
    #if kkt_err_cand / kkt_err < (1.0 - pars.kkt_reduction_factor * accept.step_size_P)
      return :success
    else
      return :not_enough_progress
    end
end


function accept_func!(accept::Class_filter_ls, iter::Class_iterate, candidate::Class_iterate, filter::Array{Class_filter,1}, pars::Class_parameters)
    status = accept_func_stable!(accept, iter, candidate, pars)
    if status == :success || status == :infeasible
      return status
    else
      return accept_func_kkt!(accept, iter, candidate, filter,  pars)
    end
end

function accept_func!(accept::Class_kkt_ls, iter::Class_iterate, candidate::Class_iterate, filter::Array{Class_filter,1}, pars::Class_parameters)
    return accept_func_kkt!(accept, iter, candidate, filter, pars)
end

function accept_func!(accept::Class_agg_ls, intial_it::Class_iterate, candidate::Class_iterate, filter::Array{Class_filter,1},  pars::Class_parameters)
    if scaled_dual_feas(candidate, pars) < get_mu(candidate) * pars.agg_protect_factor
      return :success
    else
      return :not_enough_progress
    end
end

function simple_ls(iter::Class_iterate, orginal_dir::Class_point, accept_type::Symbol, filter::Array{Class_filter,1},  pars::Class_parameters, min_step_size::Float64=0.0)

    if pars.max_step_primal_dual
      step_size_P = max_step_primal_dual(iter, orginal_dir, pars.fraction_to_boundary)
    else
      step_size_P = max_step_primal(iter, orginal_dir, pars.fraction_to_boundary)
    end

    if accept_type == :accept_stable
      accept_obj = Class_stable_ls(iter, orginal_dir, pars)
    elseif accept_type == :accept_aggressive
      accept_obj = Class_agg_ls(iter, orginal_dir, pars)
    elseif accept_type == :accept_filter
      accept_obj = Class_filter_ls(iter, orginal_dir, pars)
    elseif accept_type == :accept_kkt
      accept_obj = Class_kkt_ls(iter, orginal_dir, pars)
    else
      error("acceptance function not defined")
    end

    if accept_obj.predict_red >= 0.0
        my_warn("predicted reduction non-negative")
        #accept_obj.num_steps = 0
        #return :predict_red_non_negative, iter, accept_obj
    end

    for i = 1:pars.ls_num_backtracks
      if step_size_P >= min_step_size
        candidate, is_feas, step_size_D = move(iter, orginal_dir, step_size_P, pars)
        #is_feasible(candidate, pars.comp_feas)


        accept_obj.step_size_P = step_size_P
        accept_obj.step_size_D = step_size_D
        accept_obj.num_steps = i
        if is_feas
            status = accept_func!(accept_obj, iter, candidate, filter, pars)

            if status == :success
              return :success, candidate, accept_obj
            elseif status == :predict_red_non_negative
              return status, iter, accept_obj
            end
        end

        step_size_P *= pars.ls_backtracking_factor
      else
        return :min_α, iter, accept_obj
      end
    end

    return :max_ls_it, iter, accept_obj
end



function eigenvector_ls(iter::Class_iterate, orginal_dir::Class_point, pars::Class_parameters)
    step_size = 1.0

    best_candidate = iter;
    intial_val = eval_merit_function(iter, pars)
    best_val = intial_val

    max_it = 30;

    i = 1;
    for i = 1:max_it
      dir_pos = scale_direction(orginal_dir, step_size)
      dir_neg = scale_direction(orginal_dir, -step_size)
      candidate_pos = move(iter, dir_pos, pars)
      candidate_neg = move(iter, dir_neg, pars)

      better = false

      new_val = eval_merit_function(candidate_pos, pars)
      if new_val < best_val
        best_candidate = candidate_pos
        best_val = new_val
        better = true
      end

      new_val = eval_merit_function(candidate_neg, pars)
      if new_val < best_val
        best_candidate = candidate_neg
        best_val = new_val
        better = true
      end

      if !better
        break
      end

      step_size *= 2.0
    end

    @show step_size, norm(orginal_dir.x), best_val - intial_val

    if i == max_it
      my_warn("max it reached for eig search")
    end

    return :success, best_candidate
end
