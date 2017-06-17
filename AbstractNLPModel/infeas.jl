using NLPModels
importall NLPModels

#=function copy_meta(meta::AbstractNLPModelMeta)
    return NLPModelMeta(::Any; x0, lvar, uvar, nbv, niv, nlvb, nlvo, nlvc, nlvbi, nlvci, nlvoi, nwv, ncon, y0, lcon, ucon, nnzo, nnzj, nnzh, lin, nln, nnet, lnet, nlin, nnln, nnnet, nlnet, minimize, nlo, islp, name)
#end=#

type SimpleMeta <: AbstractNLPModelMeta
    ncon::Int64
    nvar::Int64

    x0::Array{Float64,1}
    lvar::Array{Float64,1}
    uvar::Array{Float64,1}
    lcon::Array{Float64,1}
    ucon::Array{Float64,1}

    ifix::Array{Int64,1}
    ifree::Array{Int64,1}
    ilow::Array{Int64,1}
    iupp::Array{Int64,1}

    minimize::Bool

    function SimpleMeta()
      return new()
    end
end

type Class_infeas_NLP <: AbstractNLPModel
    meta::AbstractNLPModelMeta
    inner_nlp::AbstractNLPModel
    fstar::Float64
    c::Array{Float64,1}

    function Class_infeas_NLP(nlp::AbstractNLPModel, fstar::Float64, c::Array{Float64,1})
        this = new()

        this.inner_nlp = nlp
        meta = SimpleMeta()
        meta.x0 = nlp.meta.x0
        meta.ncon = nlp.meta.ncon + 1
        meta.nvar = nlp.meta.nvar

        meta.lvar = nlp.meta.lvar
        meta.uvar = nlp.meta.uvar
        meta.lcon = [nlp.meta.lcon; 0.0]
        meta.ucon = [nlp.meta.ucon; Inf]

        meta.ifix = nlp.meta.ifix
        meta.ifree = nlp.meta.ifree
        meta.ilow = nlp.meta.ilow
        meta.iupp = nlp.meta.iupp

        meta.minimize = nlp.meta.minimize

        this.meta = meta

        this.fstar = fstar
        this.c = c

        return this
    end
end

function obj(nlp::Class_infeas_NLP, x::Array{Float64,1})
    return dot(nlp.c, x)
end

function grad(nlp::Class_infeas_NLP, x::Array{Float64,1})
    return nlp.c
end

function cons(nlp::Class_infeas_NLP, x::Array{Float64,1})
    return [cons(nlp.inner_nlp, x); nlp.fstar - obj(nlp.inner_nlp, x)]
end

function jac(nlp::Class_infeas_NLP, x::Array{Float64,1})
    return [jac(nlp.inner_nlp, x); -grad(nlp.inner_nlp, x)'];
end

function hess(nlp::Class_infeas_NLP, x::Array{Float64,1};  kwargs...)
  # obj_weight=1.0::Float64, y=false::Array{Float64,1})
    keywords = Dict(key=>value for (key, value) in kwargs)
    #keywords = Dict(zip(kwargs))
    #@show keywords
    y = keywords[:y]
    return hess(nlp.inner_nlp, x, obj_weight=y[end], y=y[1:(end-1)])
end

function hess_coord(nlp::Class_infeas_NLP, x::Array{Float64,1};  kwargs...)
    keywords = Dict(key=>value for (key, value) in kwargs)
    y = keywords[:y]
    if :obj_weight in keys(keywords)
      obj_weight = keywords[:obj_weight]
    else
      obj_weight = 1.0
    end

    return hess_coord(nlp.inner_nlp, x, obj_weight=y[end], y=y[1:(end-1)])
end

function jac_coord(nlp::Class_infeas_NLP, x::Array{Float64,1};  kwargs...)
    return jac_coord(nlp, x)
end
