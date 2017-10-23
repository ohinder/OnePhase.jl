type Class_line_search_parameters
  # fill in
end

type Class_IPM_parameters
  # fill in
end

type Class_termination_parameters

end

type Class_delta_parameters
    max::Float64
    max_it::Int64
    start::Float64
    dec::Float64
    inc::Float64
    zero::Float64 #get_mu(iter) / ((1e2 + norm(get_x(iter),Inf)) * 1e2)
    min::Float64

    function Class_delta_parameters()
        this = new()
        this.max = 10.0^(60.0)
        this.max_it = 200
        this.start = 1e-6
        this.zero = 0.0
        this.min = 1e-12
        this.inc = 8.0
        this.dec = 1.0 / pi

        return this
    end
end

abstract abstract_reduct_factors;

type Class_parameters
    output_level::Int64
    debug_mode::Int64
    throw_error_nans::Bool

    # init
    start_satisfying_bounds::Bool
    mu_primal_ratio::Float64
    init_style::Symbol

    # Class_termination_parameters
    max_it::Int64
    MAX_TIME::Float64
    tol::Float64
    tol_dual_abs::Float64
    tol_infeas::Float64

    # line search
    kkt_reduction_factor::Float64
    predict_reduction_factor::Float64
    predict_reduction_factor_MAX::Float64
    predict_reduction_eigenvector_threshold::Float64
    fraction_to_boundary_predict::Float64
    fraction_to_boundary::Float64
    fraction_to_boundary_predict_exp::Float64
    ls_backtracking_factor::Float64
    ls_num_backtracks::Int64
    ls_mode_agg::Symbol
    agg_protect_factor::Float64
    move_primal_seperate_to_dual::Bool
    dual_ls::Int64
    max_step_primal_dual::Bool
    s_update::Symbol
    stb_before_agg::Bool
    mu_update::Symbol

    # IPM GENERAL
    inertia_test::Bool
    aggressive_dual_threshold::Float64
    max_it_corrections::Int64
    dual_scale_threshold::Float64
    dual_scale_mode::Symbol
    threshold_type::Symbol
    lag_grad_test::Bool
    comp_feas::Float64
    comp_feas_agg_inf::Float64
    comp_feas_agg::Float64
    min_step_size_stable::Float64
    min_step_size_agg_ratio::Float64
    ls_mode_stable_trust::Symbol
    ls_mode_stable_delta_zero::Symbol
    ls_mode_stable_correction::Symbol
    LS_non_negative_predicted_gain::Bool
    filter_type::Symbol
    use_delta_s::Bool
    adaptive_mu::Symbol
    pause_primal::Bool
    eigen_search::Bool
    trust_region::Bool
    kkt_include_comp::Bool
    primal_bounds_dual_feas::Bool

    proximal_style::Symbol
    use_prox::Bool
    use_reg::Bool
    ##
    x_norm_penalty::Float64
    a_norm_penalty::Float64

    # SADDLE PROBLEM
    ItRefine_BigFloat::Bool
    ItRefine_Num::Int64
    saddle_err_tol::Float64
    kkt_solver_type::Symbol
    linear_solver_type::Symbol
    linear_solver_safe_mode::Bool
    move_type::Symbol
    agg_eta::Symbol

    stable_reduct_factors::abstract_reduct_factors
    aggressive_reduct_factors::abstract_reduct_factors

    delta::Class_delta_parameters

    function Class_parameters()
        this = new()

        this.delta = Class_delta_parameters()

        # init
        #this.start_satisfying_bounds = true #true #true
        this.start_satisfying_bounds = true
        this.mu_primal_ratio = 1.0 #10.0 #1.0 #1e-3
        this.init_style = :mehotra
        #this.init_style = :old_style # SOMETHING WRONG WITH THIS


        this.aggressive_dual_threshold = 1.0 #1.0 #1.0
        this.dual_scale_threshold = 100.0;
        this.threshold_type = :mu
        #this.threshold_type = :mu_primal
        #this.threshold_type = :primal
        this.lag_grad_test = true
        this.dual_scale_mode = :scaled
        #this.dual_scale_mode = :sqrt
        #this.dual_scale_mode = :exact
        #this.dual_scale_mode = :primal_dual
        this.inertia_test = false # true
        this.max_it_corrections = 2 #3 ######
        this.comp_feas_agg_inf = Inf
        this.comp_feas = 1/100.0 #1/100.0
        this.comp_feas_agg = 1/50.0 #1/50.0 #1/70.0 #1/50.0
        #this.comp_feas = 1/20.0
        #this.comp_feas_agg = 1/10.0 #1/70.0 #1/70.0 #1/50.0
        this.min_step_size_stable = 0.5^5.0
        this.min_step_size_agg_ratio = 1e-4
        this.LS_non_negative_predicted_gain = true
        this.use_delta_s = false
        this.adaptive_mu = :none
        #this.adaptive_mu = :test7 # DEFAULT
        #this.adaptive_mu = :test11
        #this.adaptive_mu = :paper
        #this.adaptive_mu = :paper2
        this.primal_bounds_dual_feas = false

        this.pause_primal = false
        this.stb_before_agg = false
        this.eigen_search = false
        this.trust_region = false

        this.output_level = 3
        this.debug_mode = 0
        this.throw_error_nans = false

        this.tol = 1e-6
        this.tol_dual_abs = 1e-6
        this.tol_infeas = 1e-12 # ????
        this.max_it = 3000;
        #this.MAX_TIME = 30.0
        this.MAX_TIME = 60.0 * 60 # 10 minutes max time

        # LINE SEARCH
        this.kkt_reduction_factor = 0.2
        this.predict_reduction_factor = 0.1 #1e-1
        this.predict_reduction_factor_MAX = 0.3
        this.predict_reduction_eigenvector_threshold = 1e-1
        this.fraction_to_boundary = 0.1
        this.fraction_to_boundary_predict = 0.2
        this.fraction_to_boundary_predict_exp = 1.5
        this.ls_backtracking_factor = 0.5
        this.ls_num_backtracks = 60;
        ls_mode = :accept_filter
        this.ls_mode_stable_trust = ls_mode
        this.ls_mode_stable_delta_zero = ls_mode
        this.ls_mode_stable_correction = ls_mode
        this.ls_mode_agg = :accept_aggressive
        this.agg_protect_factor = 10.0^6.0
        #this.agg_eta = :affine
        this.agg_eta = :mehrotra_stb
        #this.agg_eta = :mehrotra
        #this.agg_eta = :constant
        #this.protect_factor_boundary_threshold = ...
        #this.filter_type = :default
        this.filter_type = :test2
        this.kkt_include_comp = true

        this.move_type = :primal_dual
        this.move_primal_seperate_to_dual = true
        this.dual_ls = 1
        this.max_step_primal_dual = false #false
        this.s_update = :careful # :careful :loose, use careful except for experimentation
        this.mu_update = :static #:dynamic #:static #:static #:static #:dynamic :dynamic_agg

        this.saddle_err_tol = Inf
        this.ItRefine_Num = 3
        this.ItRefine_BigFloat = false
        this.use_prox = true # i.e. modify the hessian
        this.use_reg = true # i.e. modify the gradient/phi/lag
        this.proximal_style = :fixed

        this.x_norm_penalty = 1e-8 #1e-8
        this.a_norm_penalty = 1e-4

        # Don't change these parameters except for experimentation
        this.stable_reduct_factors = Reduct_stable()
        this.aggressive_reduct_factors = Reduct_affine()

        if true
          this.kkt_solver_type = :schur #_direct
          this.linear_solver_type = :julia
          #this.linear_solver_type = :mumps
        else
          this.kkt_solver_type = :symmetric
          this.linear_solver_type = :mumps
        end
        this.linear_solver_safe_mode = false #true

        return this
    end
end
