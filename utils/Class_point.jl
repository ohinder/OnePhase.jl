using Calculus

type Class_cache
    function Class_cache()
      return new()
    end
end

abstract abstract_nlp;

type Class_point
    x::Array{Float64,1}
    y::Array{Float64,1}
    s::Array{Float64,1}
    mu::Float64
    function Class_point(x::Array{Float64,1},y::Array{Float64,1},s::Array{Float64,1},mu::Float64)
        @assert(length(s) == length(y))
        return new(x, y, s, mu)
    end
end

type Class_point_v2
    x::Array{Float64,1}

    ncon::Int64

    i_l::Array{Int64,1} # variable to constraints
    y_l::Array{Float64,1}
    s_l::Array{Float64,1}

    i_u::Array{Int64,1}
    y_u::Array{Float64,1}
    s_u::Array{Float64,1}

    mu::Float64
    function Class_point_v2(nvar::Int64, ncon::Int64)
        @assert(length(s) == length(y))
        return new(x, y, s, mu)
    end
end

function get_y_divide_s(p::Class_point_v2)
    y_divide_s_l = point.y_l ./ point.s_l
    y_divide_s_u = point.y_u ./ point.s_u

    result = zeros(p.ncon)
    result[p.i_l] += y_divide_s_l
    result[p.i_u] += y_divide_s_u

    return result
end

function get_net_y(p::Class_point_v2)
    y = zeros(p.ncon)
    y[p.i_l] += p.y_l
    y[p.i_u] += p.y_u

    return y
end

function zero_point(dim::Int64,ncon::Int64)
    return Class_point(zeros(dim),zeros(ncon),zeros(ncon),0.0)
end

type Class_local_info
    delta::Float64
end

type Class_iterate
    point::Class_point
    nlp::abstract_nlp
    cache::Class_cache
    local_info::Class_local_info

    function Class_iterate(intial_point::Class_point, nlp::abstract_nlp, local_info::Class_local_info)
      return new(intial_point, nlp, Class_cache(), local_info)
    end
end

import Base.copy

function copy(it::Class_iterate)
   new_point = deepcopy(it.point)
   new_it = Class_iterate(new_point, it.nlp, it.local_info)
end

function validate(it::Class_iterate)
   approx_grad = Calculus.gradient(it.nlp.eval_f, it.point.x)
   actual_grad = it.nlp.eval_grad_lag(it.point.x, zeros(length(it.point.y)))
   @assert(norm(approx_grad -  actual_grad, 2) < 1e-3)

   actual_jac = it.nlp.eval_jac(it.point.x)
   for i = 1:length(it.point.s)
       a_i = x::Vector -> it.nlp.eval_a(x)[i]
       approx_grad_a_i = Calculus.gradient(a_i, it.point.x)
       @assert(norm(actual_jac[i,:] - approx_grad_a_i,2) < 1e-3)
   end
end

function shrink_direction!(dir::Class_point, by::Float64)
    dir.x = dir.x * by
    dir.y = dir.y * by
    dir.s = dir.s * by
    dir.mu = dir.mu * by;
end


function move(it::Class_iterate, dir::Class_point)
    new_it = copy(it)
    new_it.point.mu += dir.mu
    new_it.point.x += dir.x
    new_a = eval_a(new_it)
    new_it.point.s = new_a + (new_it.point.mu/it.point.mu) * (it.point.s - eval_a(it))
    #new_it.point.s += dir.s
    new_it.point.y += dir.y

    #new_it.point.y = new_it.point.mu ./ new_it.point.s


    return new_it;
end

function dim(it::Class_iterate)
    return length(it.point.x)
end

function ncon(it::Class_iterate)
    return length(it.point.s)
end

function get_delta(it::Class_iterate)
    return it.local_info.delta
end

function set_delta(it::Class_iterate, delta::Float64)
    it.local_info.delta = delta;
end

function get_x(it::Class_iterate)
    return it.point.x
end

function get_y(it::Class_iterate)
    return it.point.y
end

function get_s(it::Class_iterate)
    return it.point.s
end

function get_mu(it::Class_iterate)
    return it.point.mu
end
