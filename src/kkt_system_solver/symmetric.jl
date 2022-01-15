using LinearAlgebra
mutable struct Symmetric_KKT_solver <: abstract_KKT_system_solver
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
    schur_diag::Array{Float64,1}
    true_x_diag::Array{Float64,1}

    ready::Symbol #

    # Symmetric_KKT_solver only
    #true_diag::Array{Float64,1}
    Q::SparseMatrixCSC{Float64,Int64}
    #M::SparseMatrixCSC{Float64,Int64}
    #K::SparseMatrixCSC{Float64,Int64}
    #current_it::Class_iterate
    #reduct_factors::Class_reduction_factors

    function Symmetric_KKT_solver()
        ## reworks system to form of
        ## [[H A']]
        ## [[A D]].
        ##
      return new()
    end
end

function form_system!(kkt_solver::Symmetric_KKT_solver, iter::Class_iterate, timer::class_advanced_timer)
    start_advanced_timer(timer, "symmetric/form_system");
    J = get_jac(iter) #.cache.J
    J_T = get_jac_T(iter) #.cache.J_T
    #B = spdiagm(-iter.point.s ./ iter.point.y);
    B = sparse(Diagonal(-iter.point.s ./ iter.point.y))
    M = [[get_lag_hess(iter) J_T]; [J B]];

    pause_advanced_timer(timer, "symmetric/form_system");

    start_advanced_timer(timer, "symmetric/allocate");
    kkt_solver.Q = M;
    kkt_solver.factor_it = iter;
    kkt_solver.schur_diag = compute_schur_diag(iter)
    kkt_solver.true_x_diag = diag(M)[1:dim(iter)]
    kkt_solver.ready = :system_formed
    pause_advanced_timer(timer, "symmetric/allocate");

end

function factor_implementation!(kkt_solver::Symmetric_KKT_solver, timer::class_advanced_timer)
    return ls_factor!(kkt_solver.ls_solver, kkt_solver.Q, dim(kkt_solver.factor_it), ncon(kkt_solver.factor_it), timer)
end

function compute_direction_implementation!(kkt_solver::Symmetric_KKT_solver, timer::class_advanced_timer)
    start_advanced_timer(timer, "symmetric")
    y_org = get_y(kkt_solver.factor_it)
    s_org = get_s(kkt_solver.factor_it)
    rhs = kkt_solver.rhs;

    symmetric_rhs = [rhs.dual_r; rhs.primal_r + rhs.comp_r ./ y_org];
    dir_x_and_y = ls_solve(kkt_solver.ls_solver, symmetric_rhs, timer);

    shift_vector = [ones(dim(kkt_solver.factor_it));
    zeros(ncon(kkt_solver.factor_it))] * get_delta(kkt_solver.factor_it);

    dir = kkt_solver.dir;
    dir.x = dir_x_and_y[1:length(rhs.dual_r)];
    dir.y = -dir_x_and_y[(1+length(rhs.dual_r)):end];
    #dir.s = ( rhs.comp_r - dir.y .* s_org ) ./ y_org
    dir.s = eval_jac_prod(kkt_solver.factor_it, dir.x) - rhs.primal_r

    check_for_nan(dir)

    start_advanced_timer(timer, "symmetric/kkt_err");
    update_kkt_error!(kkt_solver, Inf, timer)
    pause_advanced_timer(timer, "symmetric/kkt_err");
    pause_advanced_timer(timer, "symmetric")
end

function update_delta_vecs!(kkt_solver::Symmetric_KKT_solver, delta_x_vec::Array{Float64,1}, delta_s_vec::Array{Float64,1}, timer::class_advanced_timer)
    start_advanced_timer(timer, "symmetric")
    start_advanced_timer(timer, "symmetric/delta_vecs")
    kkt_solver.delta_x_vec = delta_x_vec
    kkt_solver.delta_s_vec = delta_s_vec

    if sum(abs.(delta_s_vec)) > 0.0
        error("not implemented")
    else
        for i = 1:length(delta_x_vec)
            kkt_solver.Q[i,i] = kkt_solver.true_x_diag[i] + delta_x_vec[i]
        end
    end

    kkt_solver.ready = :delta_updated
    pause_advanced_timer(timer, "symmetric/delta_vecs")
    pause_advanced_timer(timer, "symmetric")
end
