type Class_reduction_factors
    P::Float64
    D::Float64
    mu::Float64
    function Class_reduction_factors(P::Float64,D::Float64,mu::Float64)
      return new(P,D,mu)
    end
    function Class_reduction_factors()
      return new(NaN,NaN,NaN)
    end
end

type System_rhs
    dual_r::Array{Float64,1}
    primal_r::Array{Float64,1}
    comp_r::Array{Float64,1}
    
    function System_rhs()
        return new()
    end

    function System_rhs(it::Class_iterate)
        return new(zeros(dim(it)), zeros(ncon(it)), zeros(ncon(it)))
    end

    function System_rhs(it::Class_iterate, reduct::Class_reduction_factors)
      dual_target = -eval_grad_lag(it) * (1.0 - reduct.D)
      primal_target = -eval_primal_residual(it) * (1.0 - reduct.P)
      mu_target = get_mu(it) * reduct.mu
      s = get_s(it)
      y = get_y(it)

      return new(dual_target, primal_target, mu_target - s .* y);
    end
end

import Base.LinAlg.norm
function norm(rhs::System_rhs, p::Float64)
    return norm([rhs.dual_r, rhs.primal_r, rhs.comp_r], p)
end
