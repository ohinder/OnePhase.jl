# what is the difference between schur direct and schur????

type Schur_KKT_solver_direct <: abstract_schur_solver
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
    current_it::Class_iterate
    reduct_factors::Class_reduction_factors

    function Schur_KKT_solver_direct()
      this = new()
      this.ready = :not_ready

      return this
    end
end

function compute_direction_implementation!(kkt_solver::Schur_KKT_solver_direct, timer::class_advanced_timer)
    start_advanced_timer(timer, "SCHUR")

    it = kkt_solver.current_it
    reduct = kkt_solver.reduct_factors
    pars = kkt_solver.pars

    #∇a = get_jac(it);
    y = get_y(it);
    s = get_s(it);

    rhs = kkt_solver.rhs

    #r1 + ∇a_org' * (( r3 + (r2 .* y_org) ) ./ s_org)
    start_advanced_timer(timer, "SCHUR/rhs");
    symmetric_primal_rhs = rhs.primal_r + rhs.comp_r ./ y
    Σ_vec = ( y ./ s )
    y_ = rhs.primal_r .* Σ_vec + rhs.comp_r ./ s
    schur_rhs = rhs.dual_r + eval_jac_T_prod(it, y_)
    pause_advanced_timer(timer, "SCHUR/rhs");

    dir = kkt_solver.dir
    dir.x = solver_schur_rhs(schur_rhs, kkt_solver, timer)
    dir.y = -(eval_jac_prod(it, dir.x) - symmetric_primal_rhs) .* Σ_vec
    dir.s = ( rhs.comp_r - dir.y .* s ) ./ y
    #dir.s = get_jac(it) * dir.x - rhs.primal_r
    #dir.y = ( rhs.comp_r - dir.s .* get_y(it) ) ./ get_s(it) # (mu_target - s_cur .* y_cur)

    check_for_nan(dir)

    start_advanced_timer(timer, "SCHUR/kkt_err");
    update_kkt_error!(kkt_solver, Inf, timer)
    pause_advanced_timer(timer, "SCHUR/kkt_err");
    pause_advanced_timer(timer, "SCHUR")
end
