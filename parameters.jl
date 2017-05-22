type Class_parameters
    kkt_reduction_factor::Float64
    predict_reduction_factor::Float64
    predict_reduction_factor_MAX::Float64
    predict_reduction_eigenvector_threshold::Float64
    aggressive_dual_threshold::Float64
    max_it_corrections::Int64
    dual_scale_threshold::Float64
    tol::Float64
    saddle_err_tol::Float64
    comp_feas::Float64
    comp_feas_agg::Float64
    min_step_size_stable::Float64
    min_step_size_correction::Float64
    output_level::Int64
    max_it::Int64
    MAX_TIME::Float64
    num_iterative_refinements::Int64
    fraction_to_boundary::Float64
    ls_backtracking_factor::Float64
    ls_num_backtracks::Int64
    kkt_solver_type::Symbol
    linear_solver_type::Symbol
    linear_solver_safe_mode::Bool
    move_type::Symbol
    agg_protect_factor::Float64
    ls_mode_stable::Symbol
    ls_mode_agg::Symbol
    use_delta_s::Bool
    move_primal_seperate_to_dual::Bool
    max_step_primal_dual::Bool
    s_update::Symbol

    function Class_parameters()
        this = new()
        this.kkt_reduction_factor = 0.5
        this.predict_reduction_factor = 1e-3 #1e-1
        this.predict_reduction_factor_MAX = 0.3
        this.predict_reduction_eigenvector_threshold = 1e-1
        this.aggressive_dual_threshold = 1.0
        this.dual_scale_threshold = 1.0;
        this.max_it_corrections = 1
        this.tol = 1e-6
        this.comp_feas = 1e-2
        this.comp_feas_agg = 5e-1
        this.saddle_err_tol = 0.5
        this.min_step_size_stable = 1e-2
        this.min_step_size_correction = 1e-1
        this.output_level = 3
        this.max_it = 3000;
        this.MAX_TIME = 60.0 * 10 # 10 minutes max time
        this.num_iterative_refinements = 10
        this.fraction_to_boundary = 0.05
        this.ls_backtracking_factor = 0.5
        this.ls_num_backtracks = 60;
        this.ls_mode_stable = :accept_filter
        this.ls_mode_agg = :accept_aggressive
        this.move_type = :primal_dual
        this.move_primal_seperate_to_dual = true
        this.max_step_primal_dual = false
        this.use_delta_s = false
        this.agg_protect_factor = 1e4
        this.s_update = :careful # :careful :loose, use careful except for experimentation

        if true
          this.kkt_solver_type = :schur
          this.linear_solver_type = :julia
        else
          this.kkt_solver_type = :symmetric
          this.linear_solver_type = :mumps
        end
        this.linear_solver_safe_mode = true

        return this
    end
end
