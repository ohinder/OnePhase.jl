
mutable struct Class_point
    x::Array{Float64,1}
    y::Array{Float64,1}
    s::Array{Float64,1}
    mu::Float64
    primal_scale::Float64

    function Class_point(x::Array{Float64,1},y::Array{Float64,1},s::Array{Float64,1},mu::Float64)
        #println("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", s)
        #println("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", length(s))
        #println("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", y)
        #println("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", length(y))
        @assert(length(s) == length(y))
        return new(x, y, s, mu, 1.0)
    end
end

#=function shrink_direction!(dir::Class_point, by::Float64)
    dir.x = dir.x * by
    dir.y = dir.y * by
    dir.s = dir.s * by
    dir.mu = dir.mu * by;
    dir.mu = dir.mu * by;
end=#

function zero_point(dim::Int64,ncon::Int64)
    return Class_point(zeros(dim),zeros(ncon),zeros(ncon),0.0)
end

#=
type Class_point
    x::Array{Float64,1}
    ncon::Int64

    i_l::Array{Int64,1} # variable to constraints
    y_l::Array{Float64,1}
    s_l::Array{Float64,1}

    i_u::Array{Int64,1}
    y_u::Array{Float64,1}
    s_u::Array{Float64,1}

    mu::Float64

    function Class_point(nvar::Int64, i_l::Array{Int64,1}, i_u::Array{Int64,1})
        this = new();
        this.x = zeros(nvar)
        this.i_l = i_l;
        this.i_u = i_u;

        this.y_l = zeros(length(i_l))
        this.s_l = zeros(length(i_l))

        this.s_u = zeros(length(i_u))
        this.y_u = zeros(length(i_u))

        this.mu = 0.0;

        return this
    end
end

function get_y_divide_s(p::Class_point)
    y_divide_s_l = point.y_l ./ point.s_l
    y_divide_s_u = point.y_u ./ point.s_u

    result = zeros(p.ncon)
    result[p.i_l] += y_divide_s_l
    result[p.i_u] += y_divide_s_u

    return result
end

function get_net_y(p::Class_point)
    y = zeros(p.ncon)
    y[p.i_l] += p.y_l
    y[p.i_u] += p.y_u

    return y
end

function zero_point(point::Class_point)
    new_point = deepcopy(point)
    point.x = zeros(length(point.x))

    nlb = length(point.i_l)
    point.s_l = zeros(nlb)
    point.y_l = zeros(nlb)

    nub = length(point.i_u)
    point.s_u = zeros(nub)
    point.y_u = zeros(nub)

    point.mu = 0.0

    return point
end=#
