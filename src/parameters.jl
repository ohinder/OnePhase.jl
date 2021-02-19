@compat abstract type abstract_reduct_factors end
@compat abstract type abstract_pars end

type Class_kkt_solver_options <: abstract_pars
    ItRefine_Num::Int64
    ItRefine_BigFloat::Bool
    saddle_err_tol::Float64
    kkt_solver_type::Symbol
    kkt_system_rescale::Symbol
    linear_solver_type::Symbol
    linear_solver_safe_mode::Bool
    linear_solver_recycle::Bool
    ma97_u::Float64
    stable_reduct_factors::abstract_reduct_factors
    aggressive_reduct_factors::abstract_reduct_factors

    function Class_kkt_solver_options()
      this = new()

      this.ItRefine_Num = 3
      this.ItRefine_BigFloat = false
      this.saddle_err_tol = Inf

      this.kkt_system_rescale = :none
      this.ma97_u = 1e-8
      #this.kkt_system_rescale = :u_only
      #this.kkt_system_rescale = :u_and_x
      #this.kkt_system_rescale = :u_and_x
      if true
        this.kkt_solver_type = :schur #_direct
        this.linear_solver_type = :julia
        #this.linear_solver_type = :mumps
      else
        this.kkt_solver_type = :symmetric
        this.linear_solver_type = :mumps
      end
      this.linear_solver_safe_mode = false #true
      this.linear_solver_recycle = false

      # Don't change these parameters except for experimentation
      this.stable_reduct_factors = Reduct_stable()
      this.aggressive_reduct_factors = Reduct_affine()

      return this
    end
end

type Class_line_search_parameters <: abstract_pars
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

  # to add
  comp_feas::Float64
  comp_feas_agg::Float64
  min_step_size_stable::Float64
  min_step_size_agg_ratio::Float64

  function Class_line_search_parameters()
    this = new()
    # LINE SEARCH
    this.kkt_reduction_factor = 0.2
    this.kkt_include_comp = true
    this.filter_type = :test2

    this.predict_reduction_factor = 0.1

    this.fraction_to_boundary = 0.1 #0.1
    this.fraction_to_boundary_predict = 0.2 #0.25
    this.fraction_to_boundary_linear = 0.1
    this.fraction_to_boundary_predict_exp = 0.5

    this.backtracking_factor = 0.5
    this.num_backtracks = 60;

    this.agg_gamma = :mehrotra_stb
    #this.agg_gamma = :affine

    this.move_primal_seperate_to_dual = true
    this.dual_ls = 1

    this.comp_feas = 1/100.0
    this.comp_feas_agg = 1/50.0
    this.min_step_size_stable = 0.5^5.0
    this.min_step_size_agg_ratio = 1e-4

    return this
  end
end

type Class_IPM_parameters
  # fill in
end

type Class_termination_parameters <: abstract_pars
    max_it::Int64
    max_time::Float64
    tol_opt::Float64
    tol_unbounded::Float64
    grad_max::Float64
    tol_inf_1::Float64
    tol_inf_2::Float64
    dual_scale_threshold::Float64
    dual_scale_mode::Symbol

    function Class_termination_parameters()
        this = new()

        this.max_it = 3000
        this.max_time = 60.0^2
        this.tol_opt = 1e-6
        this.tol_unbounded = 1e-12
        this.grad_max = Inf
        this.tol_inf_1 = 1e-3
        this.tol_inf_2 = 1e-6
        this.dual_scale_threshold = 100.0;
        this.dual_scale_mode = :max_dual

        return this
    end
end

type Class_delta_parameters <: abstract_pars
    max::Float64
    #max_it::Int64
    start::Float64
    dec::Float64
    inc::Float64
    zero::Float64 #get_mu(iter) / ((1e2 + norm(get_x(iter),Inf)) * 1e2)
    min::Float64

    function Class_delta_parameters()
        this = new()
        this.max = 10.0^(50.0)
        #this.max_it = 200
        this.start = 1e-6
        this.zero = 0.0
        this.min = 1e-12
        this.inc = 8.0
        this.dec = 1.0 / pi

        return this
    end
end


type Class_init_parameters <: abstract_pars
    mu_scale::Float64
    mehotra_scaling::Bool
    init_style::Symbol
    start_satisfying_bounds::Bool
    dual_threshold::Float64
    linear_scale::Float64
    nl_ineq_scale::Float64
    nl_eq_scale::Float64
    dual_max::Float64
    dual_min::Float64

    function Class_init_parameters()
        this = new()
        #mode = :experimentation
        mode = :standard
        if mode == :standard
          this.mu_scale = 1.0
          this.mehotra_scaling = true
          this.init_style = :gertz # #
          this.dual_threshold = 1.0
          this.start_satisfying_bounds = true
          this.linear_scale = 1.0
          this.nl_ineq_scale = 1.0
          this.nl_eq_scale = 1.0
          this.dual_max = 1e3
          this.dual_min = 1e-2
        else
          this.mu_scale = 1.0
          this.mehotra_scaling = false
          this.init_style = :mehrotra
          this.dual_threshold = 1.0
          this.start_satisfying_bounds = true
          this.linear_scale = 1e-2
          this.nl_ineq_scale = 1e-1
          this.nl_eq_scale = 1.0
        end

        return this
    end
end

type Class_testing <: abstract_pars
    response_to_failure::Symbol
    function Class_testing()
        this = new()
        #this.response_to_failure = :default
        this.response_to_failure = :lag_delta_inc

        return this
    end
end

type Class_parameters <: abstract_pars
    term::Class_termination_parameters
    init::Class_init_parameters
    delta::Class_delta_parameters
    ls::Class_line_search_parameters
    kkt::Class_kkt_solver_options
    test::Class_testing


    # debugging
    output_level::Int64
    debug_mode::Int64
    throw_error_nans::Bool

    # IPM GENERAL
    aggressive_dual_threshold::Float64
    max_it_corrections::Int64

    superlinear_theory_mode::Bool
    agg_protection_factor::Float64

    kkt_include_comp::Bool
    primal_bounds_dual_feas::Bool
    a_norm_penalty::Float64

    eps_mach::Float64

    function Class_parameters()
        this = new()

        ############################
        ### groups of parameters ###
        ############################
        this.delta = Class_delta_parameters()
        this.term = Class_termination_parameters()
        this.init = Class_init_parameters()
        this.ls = Class_line_search_parameters()
        this.kkt = Class_kkt_solver_options()
        this.test = Class_testing()

        ########################
        ## general parameters ##
        ########################

        # debugging
        this.output_level = 2
        this.debug_mode = 0
        this.throw_error_nans = false

        # switching condition
        this.aggressive_dual_threshold = 1.0 # kappa_{1}
        this.primal_bounds_dual_feas = false

        # algorithm parameters
        this.max_it_corrections = 2
        this.superlinear_theory_mode = false
        this.agg_protection_factor = 0.9

        # penalty parameter size, changes log barrier to stop iterates growing too big when feasible region is unbounded
        # e.g. min x s.t. x, y >= 0, without any changes y -> \infty
        this.a_norm_penalty = 1e-4

        # machine precision
        this.eps_mach = 1e-16

        return this
    end
end

#::IOStream
function write_pars(stream, par::Class_parameters)
    write(stream, "PAR \t\t\t\t\t\t VALUE \n")

    for fieldname in fieldnames(par)
        fieldval = getfield(par, fieldname)
        if isa(fieldval,abstract_pars)
          write(stream, "$(pd(fieldname,40)) \n")
          for subfieldname in fieldnames(fieldval)
            subfieldval = getfield(fieldval, subfieldname)
            write(stream, "\t $(pd(subfieldname,40)) \t $(subfieldval) \n")
          end
        else
          write(stream, "$(pd(fieldname,40)) \t $(fieldval) \n")
        end
    end
end
