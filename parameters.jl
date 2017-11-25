type Class_line_search_parameters
  kkt_reduction_factor::Float64
  kkt_include_comp::Bool
  filter_type::Symbol

  predict_reduction_factor::Float64

  fraction_to_boundary_predict::Float64
  fraction_to_boundary::Float64
  fraction_to_boundary_predict_exp::Float64
  fraction_to_boundary_linear::Float64

  backtracking_factor::Float64
  num_backtracks::Int64

  agg_gamma::Symbol

  move_primal_seperate_to_dual::Bool
  dual_ls::Int64

  function Class_line_search_parameters()
    this = new()
    # LINE SEARCH
    this.kkt_reduction_factor = 0.2
    this.kkt_include_comp = true
    this.filter_type = :test2

    this.predict_reduction_factor = 0.1

    this.fraction_to_boundary = 0.1 #0.1
    this.fraction_to_boundary_predict = 0.25 #0.2
    this.fraction_to_boundary_linear = 0.1
    this.fraction_to_boundary_predict_exp = 1.5

    this.backtracking_factor = 0.5
    this.num_backtracks = 60;

    this.agg_gamma = :mehrotra_stb

    this.move_primal_seperate_to_dual = true
    this.dual_ls = 1

    return this
  end
end

type Class_IPM_parameters
  # fill in
end

type Class_termination_parameters
    max_it::Int64
    max_time::Float64
    tol_opt::Float64
    tol_unbounded::Float64
    tol_inf_1::Float64
    tol_inf_2::Float64
    max_gradient::Float64

    function Class_termination_parameters()
        this = new()

        this.max_it = 3000
        this.max_time = 60.0^2
        this.tol_opt = 1e-6
        this.tol_unbounded = 1e-10
        this.tol_inf_1 = 1e-3
        this.tol_inf_2 = 1e-6
        this.max_gradient = 1e15

        return this
    end
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


type Class_init_parameters
    mu_scale::Float64
    mehotra_scaling::Bool
    init_style::Symbol
    start_satisfying_bounds::Bool
    dual_threshold::Float64
    linear_scale::Float64
    nl_scale::Float64

    function Class_init_parameters()
        this = new()
        #mode = :experimentation
        mode = :standard
        if mode == :standard
          this.mu_scale = 1.0
          this.mehotra_scaling = true
          this.init_style = :mehotra
          this.dual_threshold = 1.0
          this.start_satisfying_bounds = true
          this.linear_scale = 1.0
          this.nl_scale = 1.0
        else
          this.mu_scale = 1500
          this.mehotra_scaling = false
          this.init_style = :mehotra
          this.dual_threshold = 1.0
          this.start_satisfying_bounds = true
          this.linear_scale = 1.0
          this.nl_scale = 10.0
        end

        return this
    end
end

abstract abstract_reduct_factors;

type Class_parameters
    term::Class_termination_parameters
    init::Class_init_parameters
    delta::Class_delta_parameters
    ls::Class_line_search_parameters

    # debugging
    output_level::Int64
    debug_mode::Int64
    throw_error_nans::Bool

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
    superlinear_theory_mode::Bool

    use_delta_s::Bool
    pause_primal::Bool
    eigen_search::Bool
    trust_region::Bool
    kkt_include_comp::Bool
    primal_bounds_dual_feas::Bool
    #Fmu_scale::Float64

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
    #move_type::Symbol

    stable_reduct_factors::abstract_reduct_factors
    aggressive_reduct_factors::abstract_reduct_factors

    function Class_parameters()
        this = new()

        this.delta = Class_delta_parameters()
        this.term = Class_termination_parameters()
        this.init = Class_init_parameters()
        this.ls = Class_line_search_parameters()

        # init
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
        this.max_it_corrections = 3 #3 ######
        this.comp_feas_agg_inf = Inf
        this.comp_feas = 1/100.0 #1/100.0
        this.comp_feas_agg = 1/50.0 #1/50.0 #1/70.0 #1/50.0
        #this.comp_feas = 1/20.0
        #this.comp_feas_agg = 1/10.0 #1/70.0 #1/70.0 #1/50.0
        this.min_step_size_stable = 0.5^5.0
        this.min_step_size_agg_ratio = 1e-4
        this.use_delta_s = false
        this.primal_bounds_dual_feas = false
        #this.mu_scale = 1.0

        this.superlinear_theory_mode = false

        this.pause_primal = false
        this.eigen_search = false
        this.trust_region = false

        this.output_level = 3
        this.debug_mode = 0
        this.throw_error_nans = false

        this.saddle_err_tol = Inf
        this.ItRefine_Num = 3
        this.ItRefine_BigFloat = false
        this.use_prox = true # i.e. modify the hessian
        this.use_reg = true # i.e. modify the gradient/phi/lag
        this.proximal_style = :fixed

        this.x_norm_penalty = 0.0 #1e-8
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
