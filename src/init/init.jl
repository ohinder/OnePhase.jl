include("correct-guess.jl")
include("guess-vars.jl")
include("other.jl")
include("primal-project.jl")
include("gertz_init.jl")
include("lp.jl")

function linear_cons(nlp::Class_CUTEst, x0::Vector)
    a1 = eval_a(nlp, x0)
    d = randn(length(x0))
    a2 = eval_a(nlp, x0 + d)
    a3 = eval_a(nlp, x0 + 2 * d)

    return abs(a1 - 2 * a2 + a3) .< 1e-6
end

function mehrotra_init(nlp::Class_CUTEst, pars::Class_parameters, timer::class_advanced_timer)
    start_advanced_timer(timer, "INIT/x")
    x0 = suggested_starting_point(nlp)
    x0 += rand(length(x0)) * 1.0

    if pars.init.start_satisfying_bounds
      x = projection_onto_bounds_ipopt_style( nlp, pars, x0 )
    end
    pause_advanced_timer(timer, "INIT/x")

    a, J, g = eval_init(nlp, pars, timer, x)


    s, y = mehortra_guess( nlp, pars, timer, x, a, J, g )

    @assert(all(s .>= 0.0))

    li = linear_cons(nlp)
    nl = li .!= true
    ineq = ineq_cons(nlp)
    eq = ineq .!= true
    nl_eq = ineq .& nl
    nl_ineq = eq .& nl

    if pars.init.mehotra_scaling
      mu = Statistics.mean(s .* y)
      conWeight = ((s - a) / mu)
      #v = conWeight[conWeight .> 0.0]
      #conWeight[conWeight .> 0.0] = min.(10.0,max.(v,1e-4))
    else
      mu = (1e-6 + LinearAlgebra.norm(s,Inf) + LinearAlgebra.norm(g,Inf))
      #mu = 1e-6 + LinearAlgebra.norm(g,1) / length(s)
      conWeight = zeros(length(s))
      ais = cons_indicies(nlp)
      conWeight[ais] = 1.0

      #li = linear_cons(nlp, x)
    end

    #@show li
    #@show nl

    #@show li
    if true
      conWeight[nl_eq] *= pars.init.nl_eq_scale
      conWeight[nl_ineq] *= pars.init.nl_ineq_scale
      #@show li, nl_ineq
      conWeight[li] *= pars.init.linear_scale
      #mu *= pars.init.mu_scale
    end
    #*= pars.init.linear_scale

    @assert(all(conWeight .>= 0.0))

    iter_init = correct_guess3( nlp, pars, timer, x, a, J, g, y, mu, conWeight )
    #return correct_guess1( nlp, pars, timer, x0, a, J, g, s, y)
    #
    #@show li

    #@show li
    iter_init.frac_bd_predict[li] = pars.ls.fraction_to_boundary_linear
    iter_init.frac_bd[li] = pars.ls.fraction_to_boundary_linear

    @assert(all(iter_init.point.s .>= 0.0))
    @assert is_feasible(iter_init,pars.ls.comp_feas)

    return iter_init
end

function simple_init()

end
