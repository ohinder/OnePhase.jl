type Schur_KKT_solver <: abstract_KKT_system_solver
    # abstract_KKT_system_solver
    ls_solver::abstract_linear_system_solver
    factor_it::Class_iterate
    delta_x_vec::Array{Float64,1}
    delta_s_vec::Array{Float64,1}
    rhs::System_rhs
    dir::Class_point
    kkt_err_norm::Class_kkt_error
    rhs_norm::Float64
    pars::Class_parameters

    ready::Symbol #

    # Schur_KKT_solver only
    true_diag::Array{Float64,1}
    M::SparseMatrixCSC{Float64,Int64}
    K::SparseMatrixCSC{Float64,Int64}

    function Schur_KKT_solver()
      this = new()
      this.ready = :not_ready

      return this
    end
end


function form_system!(kkt_solver::Schur_KKT_solver, iter::Class_iterate, timer::class_advanced_timer)
    start_advanced_timer(timer, "SCHUR")
    start_advanced_timer(timer, "SCHUR/form_system");
    x = iter.point.x;
    y = iter.point.y
    s = iter.point.s

    ∇a = get_jac(iter)

    #scaling = zeros(dim(iter))
    #for i = 1:dim(iter)
    #    scaling[i] = norm(∇a[:,i], Inf)
    #end

    #@show size(∇a)
    #@show maximum(scaling), minimum(scaling)
    #tic()

    Σ = spdiagm(iter.point.y ./ iter.point.s)

    kkt_solver.M = (∇a' * Σ * ∇a) + get_lag_hess(iter);
    kkt_solver.true_diag = diag(kkt_solver.M)
    kkt_solver.factor_it = iter;
    kkt_solver.ready = :system_formed

    pause_advanced_timer(timer, "SCHUR/form_system");
    pause_advanced_timer(timer, "SCHUR");
end



function update_delta_vecs!(kkt_solver::Schur_KKT_solver, delta_x_vec::Array{Float64,1}, delta_s_vec::Array{Float64,1}, timer::class_advanced_timer)
    start_advanced_timer(timer, "SCHUR")
    start_advanced_timer(timer, "SCHUR/delta_vecs")
    kkt_solver.delta_x_vec = delta_x_vec
    kkt_solver.delta_s_vec = delta_s_vec
    ∇a = ∇a = get_jac(kkt_solver.factor_it)

    if sum(abs(delta_s_vec)) > 0.0
        kkt_solver.K = kkt_solver.M + (∇a' * spdiagm(delta_s_vec) * ∇a) + spdiagm(delta_x_vec)
    else
        kkt_solver.K = kkt_solver.M + spdiagm(delta_x_vec)
    end

    kkt_solver.ready = :delta_updated
    pause_advanced_timer(timer, "SCHUR/delta_vecs")
    pause_advanced_timer(timer, "SCHUR")
end

function factor_implementation!(kkt_solver::Schur_KKT_solver, timer::class_advanced_timer)
    return ls_factor!(kkt_solver.ls_solver, kkt_solver.K, dim(kkt_solver.factor_it), 0, timer)
end

function compute_direction_implementation!(kkt_solver::Schur_KKT_solver, timer::class_advanced_timer)
    start_advanced_timer(timer, "SCHUR")

    factor_it = kkt_solver.factor_it
    ∇a_org = get_jac(factor_it);
    y_org = get_y(factor_it);
    s_org = get_s(factor_it);
    pars = kkt_solver.pars

    rhs = kkt_solver.rhs

    #r1 + ∇a_org' * (( r3 + (r2 .* y_org) ) ./ s_org)
    start_advanced_timer(timer, "SCHUR/rhs");
    symmetric_primal_rhs = rhs.primal_r + rhs.comp_r ./ y_org
    Σ_vec = ( y_org ./ s_org )
    schur_rhs = rhs.dual_r + ∇a_org' * (rhs.primal_r .* Σ_vec + rhs.comp_r ./ s_org )
    pause_advanced_timer(timer, "SCHUR/rhs");

    #(symmetric_primal_rhs .* Σ_vec);Σ
    dir = kkt_solver.dir;

    # generalize!!!
    output_level = pars.output_level
    res_old = schur_rhs
    if output_level >= 4
      println("res", 0, " ", rd(norm(res_old,2)))
    end

    dir_x = zeros(length(dir.x));
    if pars.ItRefine_BigFloat
      dir_x = convert(Array{BigFloat,1}, dir_x)
    end

    for i = 1:pars.ItRefine_Num
        start_advanced_timer(timer, "SCHUR/iterative_refinement");
        dir_x += ls_solve(kkt_solver.ls_solver, res_old, timer)[:];

        start_advanced_timer(timer, "SCHUR/iterative_refinement/residual");
        jac_res = ∇a_org' * ( Σ_vec  .* (∇a_org * dir_x) )
        hess_res = hess_product(factor_it, dir_x) + kkt_solver.delta_x_vec .* dir_x
        res = schur_rhs - ( jac_res + hess_res )
        pause_advanced_timer(timer, "SCHUR/iterative_refinement/residual");

        if output_level >= 4
          println("res", i, " ", rd(Float64(norm(res,2))))
        end

        res_old = res
        pause_advanced_timer(timer, "SCHUR/iterative_refinement");
    end
    dir.x = dir_x

    # there are two ways to update s and y
    if true
      dir.y = -(∇a_org * dir.x - symmetric_primal_rhs) .* Σ_vec
      dir.s = ( rhs.comp_r - dir.y .* s_org ) ./ y_org
    else
      dir.s = ∇a_org * dir.x - rhs.primal_r
      dir.y = ( rhs.comp_r - dir.s .* y_org ) ./ s_org # (mu_target - s_cur .* y_cur)
    end

    check_for_nan(dir)

    start_advanced_timer(timer, "SCHUR/kkt_err");
    update_kkt_error!(kkt_solver, Inf, timer)
    pause_advanced_timer(timer, "SCHUR/kkt_err");
    pause_advanced_timer(timer, "SCHUR")
end
