function one_phase(intial_it::Class_iterate, par::Class_parameters)
    iter = intial_it;
    initialize!(par.kkt_solver)

    status = :success
    step_size_stabilization = NaN

    progress = Array{Dict{String,Any},1}()

    head_progress()
    record_progress!(0,"ini", iter, zero_point(dim(iter),ncon(iter)), NaN, NaN, progress, par)

    t = 0
    for p = 1:500
        t += 1

        stable_iter, status, step_size_stabilization, dir_stable, fact_M = stabilization_step(iter, par)

        @assert(status == :success)
        # DISPLAY
        record_progress!(t,"sta", stable_iter, dir_stable, 0.0, step_size_stabilization, progress, par)

        if scaled_dual_feas(stable_iter, par) < get_mu(stable_iter)
            if false
              t += 1
              M = form_mat(iter)
              fact_M = fact_mat(M, iter)
            end

            new_iter = stable_iter
            for i = 1:100
                new_iter, status, step_size_affine, dir_affine, fact_M = affine_step(new_iter, iter, par, fact_M, progress, par)

                aff_iter = new_iter
                # DISPLAY
                record_progress!(t, "aff", new_iter, dir_affine, 1.0, step_size_affine)

                if step_size_affine < 0.5 && i > 2 || terminate_optimal(new_iter, par)
                    break
                end
            end

            iter = new_iter
        else
            iter = stable_iter;
        end

        alg_status = terminate(iter, par)
        if alg_status != false
            return iter, alg_status
        end

        set_delta(iter, get_delta(iter) / 10.0)
    end

    return iter, :max_it
end


function stabilization_step(iter::Class_iterate, par::Class_parameters)
    M = form_mat(iter)

    for delta_it = 1:10
        fact_M = fact_mat(M, iter)

        dir_stabilization, dir_status = compute_direction(iter, iter, fact_M, par, :stabilization)

        if dir_status == :success
          new_iterate, ls_status, step_size = stabilization_line_search(iter, dir_stabilization, par)
          if ls_status == :success
              return new_iterate, :success, step_size, dir_stabilization, fact_M
          end
        end

        iter.local_info.delta = max(iter.local_info.delta * 10.0, min(iter.point.mu, 1e-3))
    end

    return iter, :failure, NaN, NaN, NaN
end

function affine_step(cur_iter, org_iter, par, fact_M)
    dir_affine, dir_status = compute_direction(cur_iter, org_iter, fact_M, par, :affine)

    if dir_status == :success
      new_iterate, ls_status, step_size = affine_line_search(cur_iter, dir_affine, par)
      if ls_status == :success
          return new_iterate, :success, step_size, dir_affine, fact_M
      end
    end

    return cur_iter, :failure, NaN, NaN, NaN
end

function stabilization_line_search(iter::Class_iterate, dir_stable::Class_point, par::Class_parameters)
    step_size = 1.0
    for i = 1:30
        candidate = move(iter, dir_stable)
        status = acceptable(iter, candidate, dir_stable, par)

        # TO DO!!!
        # make sure there is sufficient progress on complementarity!
        #&& norm(comp(candidate),Inf) < max(iter.point.mu * 0.5, norm(comp(iter),Inf))

        if status == :success
            return candidate, :success, step_size;
        end

        shrink_by = 0.5
        step_size *= shrink_by
        shrink_direction!(dir_stable, shrink_by)
    end

    return iter, :failure, NaN;
end



function affine_line_search(iter::Class_iterate, org_dir::Class_point, par::Class_parameters)
    step_size = 0.5
    dir = scale_direction(org_dir, step_size)
    intial_candidate = move(iter, dir)

    #comp_dont_increase = norm(comp(candidate),Inf) / candidate.point.mu <= (0.9 + norm(comp(iter),Inf)) / (iter.point.mu * 2)

    if is_feasible(intial_candidate, par.comp_feas)# && comp_dont_increase
        residual_step_size = 1.0 - step_size
        candidate = intial_candidate

        for i = 1:20
            new_step_size = 1.0 - residual_step_size
            dir = scale_direction(org_dir, new_step_size)
            new_candidate = move(iter, dir)

            dual_reduction = norm(eval_grad_lag(new_candidate),Inf)  / (par.tol + norm(eval_grad_lag(iter),Inf))
            mu_reduction = get_mu(new_candidate) / get_mu(iter)

            @show i, dual_reduction, mu_reduction

            if !is_feasible(new_candidate, par.comp_feas) || dual_reduction / 10.0 > mu_reduction
                return candidate, :success, step_size
            end

            step_size = new_step_size
            candidate = new_candidate

            residual_step_size *= 0.3
        end

        error("this should not occur")
    else
        for i = 1:50
            dir = scale_direction(org_dir, step_size)
            candidate = move(iter, dir)

            if is_feasible(candidate, par.comp_feas)
                return candidate, :success, step_size
            end

            step_size *= 0.5
        end
    end


    return iter, :failure, NaN;
end

function form_mat(it::Class_iterate)
    x = it.point.x;
    y = it.point.y
    s = it.point.s

    Σ = spdiagm(it.point.y ./ it.point.s)
    ∇a = eval_jac(it)
    M =  ∇a' * Σ * ∇a + eval_lag_hess(it)

    return M
end


function fact_mat(M, iter::Class_iterate)
    for i = 1:100
      try
        factorize_this = Hermitian(M + get_delta(iter) * speye(dim(iter)))
        return cholfact(factorize_this)
      catch(e)
        iter.local_info.delta = max(get_delta(iter) * 10.0, min(get_mu(iter),1e-3 ))
      end
    end

    error("Delta gets too big!")
end

function compute_direction(cur_it::Class_iterate, orginal_it::Class_iterate, factorized_mat, par::Class_parameters, direction_type::Symbol)
    s_org = orginal_it.point.s
    y_org = orginal_it.point.y

    mu_cur = cur_it.point.mu
    s_cur = cur_it.point.s
    y_cur = cur_it.point.y
    ∇a_org = eval_jac(orginal_it)



    if direction_type == :stabilization
        eta = 0.0
    elseif direction_type == :affine
        eta = 1.0
    else
      error("direction type should be stabilization or affine!")
    end


    primal_target = eval_primal_residual(cur_it) * eta;
    mu_target = mu_cur * (1.0 - eta)

    r1 = eval_grad_lag(cur_it)
    r2 = primal_target
    r3 = s_cur .* y_cur - mu_target

    dir = zero_point( dim(cur_it), ncon(cur_it) )
    dir.mu = mu_target - mu_cur

    for i = 1:6
      r1 = predicted_lag(cur_it, orginal_it, dir)
      r2 = primal_target + ∇a_org * dir.x - dir.s
      r3 = mu_target - (s_cur .* y_cur + dir.s .* y_org + dir.y .* s_org)

      rhs = r1 + ∇a_org' * (( r3 + (r2 .* y_org) ) ./ s_org)
      dir.x += (factorized_mat \ -rhs)[:]

      dir.s = ∇a_org * dir.x + primal_target
      dir.y = ( (mu_target - s_cur .* y_cur) - dir.s .* y_org ) ./ s_org


      #@show norm(predicted_lag(cur_it, orginal_it, dir),2)
    end



    if norm(predicted_lag(cur_it, orginal_it, dir),Inf) > get_mu(cur_it) * 1e-4
        @show norm(predicted_lag(cur_it, orginal_it, dir),Inf)
        #=dir2 = compute_direction_LP(cur_it, orginal_it, par)
        @show norm(predicted_lag(cur_it, orginal_it, dir2),Inf)
        @show norm(mu_target - (s_cur .* y_cur + dir2.s .* y_org + dir2.y .* s_org),2)
        @show norm(dir2.x - dir.x,2)
        @show norm(dir2.s - dir.s,2)
        @show norm(dir2.y - dir.y,2)
        #@show s_cur .* y_cur + s_dir .* y_org + y_dir .* s_org
        #@show s_cur .* y_cur + dir2.s .* y_org + dir2.y .* s_org=#
    end

    return dir, :success
end


function compute_direction_LP(cur_it::Class_iterate, orginal_it::Class_iterate, eta::Float64, par::Class_parameters)
    ∇a_org = eval_jac(orginal_it)
    n = dim(cur_it)
    m = ncon(cur_it)
    s = cur_it.point.s
    y = cur_it.point.y
    S = diagm(orginal_it.point.s)
    Y = diagm(orginal_it.point.y)

    KKT_mat =
    [[eval_lag_hess(orginal_it)  -∇a_org'         zeros(n,m)]
    [-∇a_org          zeros(m,m)  eye(m,m)]
    [zeros(m,n)   S           Y]];

    mu_target = cur_it.point.mu * (1.0 - eta)
    primal_target = eval_primal_residual(cur_it) * eta

    r1 = eval_grad_lag(cur_it)
    r2 = primal_target
    r3 = mu_target - s .* y;

    rhs = [r1; r2; r3] #- it.point.mu * ones(m)

    dir = LUfact(KKT_mat) \ -rhs

    return Class_point(dir[1:n], dir[(n+1):(n+m)], dir[(n + m + 1):(n+ 2 * m)], -eta * cur_it.point.mu)
end
