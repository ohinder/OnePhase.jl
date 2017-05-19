type Symmetric_KKT_solver <: abstract_KKT_system_solver
    ls_solver::abstract_linear_system_solver
    M::SparseMatrixCSC{Float64,Int64}
    factor_it::Class_iterate
    rhs::System_rhs
    dir::Class_point
    kkt_err_norm::Class_kkt_error
    rhs_norm::Float64
    delta::Float64
    pars::Class_parameters

    function Symmetric_KKT_solver()
      return new()
    end
end

function form_system!(kkt_solver::Symmetric_KKT_solver, iter::Class_iterate)
    start_advanced_timer("symmetric/form_system");
    A = eval_jac(iter)
    B = spdiagm(-iter.point.s ./ iter.point.y);
    M = [[eval_lag_hess(iter) A']; [A B]];

    pause_advanced_timer("symmetric/form_system");

    start_advanced_timer("symmetric/allocate");
    kkt_solver.M = M;
    kkt_solver.factor_it = iter;
    pause_advanced_timer("symmetric/allocate");

end

function factor!(kkt_solver::Symmetric_KKT_solver, shift::Float64)
    kkt_solver.delta = shift
    shift_vector = [ones(dim(kkt_solver.factor_it)); zeros(ncon(kkt_solver.factor_it))] * shift;

    return ls_factor!(kkt_solver.ls_solver, kkt_solver.M + spdiagm(shift_vector), dim(kkt_solver.factor_it), ncon(kkt_solver.factor_it))
end

function compute_direction!(kkt_solver::Symmetric_KKT_solver)
    y_org = get_y(kkt_solver.factor_it)
    s_org = get_s(kkt_solver.factor_it)
    ∇a_org = eval_jac(kkt_solver.factor_it)
    rhs = kkt_solver.rhs;

    symmetric_rhs = [rhs.dual_r; rhs.primal_r + rhs.comp_r ./ y_org];
    dir_x_and_y = ls_solve(kkt_solver.ls_solver, symmetric_rhs)[:];

    shift_vector = [ones(dim(kkt_solver.factor_it));
     zeros(ncon(kkt_solver.factor_it))] * get_delta(kkt_solver.factor_it);


    dir = kkt_solver.dir;
    dir.x = dir_x_and_y[1:length(rhs.dual_r)];
    dir.y = -dir_x_and_y[(1+length(rhs.dual_r)):end];
    #dir.s = ( rhs.comp_r - dir.y .* s_org ) ./ y_org
    dir.s = ∇a_org * dir.x - rhs.primal_r

    update_kkt_error!(kkt_solver, Inf)
end
