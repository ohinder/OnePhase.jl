type Class_parameters
    predict_reduction_factor::Float64
    tol::Float64
    comp_feas::Float64
    output_level::Int64
    max_it::Int64
    fraction_to_boundary::Float64
    ls_backtracking_factor::Float64
    ls_num_backtracks::Int64
    kkt_solver_type::Symbol
    linear_solver_type::Symbol
    linear_solver_safe_mode::Bool

    function Class_parameters()
        this = new()
        this.predict_reduction_factor = 1e-3
        this.tol = 1e-6
        this.comp_feas = 1e-2
        this.output_level = 3
        this.max_it = 2000;
        this.fraction_to_boundary = 0.1
        this.ls_backtracking_factor = 0.5
        this.ls_num_backtracks = 60;
        this.kkt_solver_type = :schur
        this.linear_solver_type = :julia
        this.linear_solver_safe_mode = true

        return this
    end
end
