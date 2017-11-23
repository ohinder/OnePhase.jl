include("correct-guess.jl")
include("guess-vars.jl")
include("other.jl")
include("primal-project.jl")

function linear_cons(nlp::Class_CUTEst, x0::Vector)
    a1 = eval_a(nlp, x0)
    d = randn(length(x0))
    a2 = eval_a(nlp, x0 + d)
    a3 = eval_a(nlp, x0 + 2 * d)

    return abs(a1 - 2 * a2 + a3) .< 1e-6
end

function init(nlp::Class_CUTEst, pars::Class_parameters, timer::class_advanced_timer)
    start_advanced_timer(timer, "INIT/x")
    x0 = suggested_starting_point(nlp)
    pause_advanced_timer(timer, "INIT/x")

    if pars.init.start_satisfying_bounds
      start_advanced_timer(timer, "INIT/projection_onto_bounds")
      x = projection_onto_bounds_ipopt_style( nlp, pars, x0 )
      pause_advanced_timer(timer, "INIT/projection_onto_bounds")
    end

    a, J, g = eval_init(nlp, pars, timer, x)


    s, y = mehortra_guess( nlp, pars, timer, x, a, J, g )

    if pars.init.mehotra_scaling
      mu = mean(s .* y)
      conWeight = ((s - a) / mu) / pars.init.mu_scale
    else
      mu = 1e-6 + norm(g,1) / length(s)
      conWeight = zeros(length(s))
      ais = cons_indicies(nlp)
      conWeight[ais] = 1.0 / pars.init.mu_scale
      #li = linear_cons(nlp, x)
      #conWeight[li] /= 10.0

    end

    @assert(all(conWeight .>= 0.0))

    iter_init = correct_guess2( nlp, pars, timer, x, a, J, g, y, mu, conWeight )
    #return correct_guess1( nlp, pars, timer, x0, a, J, g, s, y)
    #
    #@show li
    li = linear_cons(nlp)
    #@show li
    iter_init.frac_bd_predict[li] = pars.ls.fraction_to_boundary_linear
    iter_init.frac_bd[li] = pars.ls.fraction_to_boundary_linear
    return iter_init
end
