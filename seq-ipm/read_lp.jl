using MAT, JuMP

function read_lp(name::String)
    lp_data = matopen("../netlib/$(name).mat")

    A = read(lp_data,"A")
    b = read(lp_data,"b")[:]
    c = read(lp_data,"c")[:]

    lbounds = nothing
    ubounds = nothing
    try
      lbounds = read(lp_data,"lbounds")[:]
      ubounds = read(lp_data,"ubounds")[:]
    catch(e)
      lbounds = zeros(length(c))
      ubounds = Inf * ones(length(c))
    end

    return A, b, c, lbounds, ubounds
end

function read_lp_into_Class_QP(name::String)
  A, b, c, lbounds, ubounds = read_lp(name)

  threshold = 1e20
  li = (lbounds .> -threshold)
  ui = (ubounds .< threshold)

  n = length(c)
  Ahat = [A; -A; speye(n)[li,:]; -speye(n)[ui,:]];
  bhat = [b; -b; lbounds[li]; -ubounds[ui]];

  return Class_QP(bhat, Ahat, c)
end
function read_lp_into_JuMP(name::String, perturb::Float64=0.0, pretend_is_nlp::Bool=false)
  A, b, c, lbounds, ubounds = read_lp(name)

  read_lp_into_JuMP(A, b, c, lbounds, ubounds, perturb, pretend_is_nlp)
end

function read_lp_into_JuMP(A, b, c, lbounds, ubounds, perturb::Float64=0.0, pretend_is_nlp::Bool=false)
  threshold = 1e20
  li = (lbounds .> -threshold)
  ui = (ubounds .< threshold)

  mod = Model()

  n = length(c)
  x = Array{JuMP.Variable,1}(n)
  for i = 1:n
    if true == li[i] && true == ui[i]
      x[i] = @variable(mod, lowerbound = lbounds[i] - perturb, upperbound = ubounds[i] + perturb)
    elseif true == li[i]
      x[i] = @variable(mod, lowerbound = lbounds[i] - perturb)
    elseif true == ui[i]
      x[i] = @variable(mod, upperbound = ubounds[i] + perturb)
    else
      x[i] = @variable(mod)
    end
  end

  if pretend_is_nlp
    @NLobjective(mod, Min, sum(c[i] * x[i] for i = 1:n))
  else
    @objective(mod, Min, sum(c[i] * x[i] for i = 1:n))
  end
  if length(b) > 0
    @constraint(mod, A * x .== b)
  end

  return mod
end

function run_IPOPT(A, b, c, lbounds, ubounds, max_iter, bound_relax_factor,mehrotra_algorithm)
    m = read_lp_into_JuMP(A, b, c, lbounds, ubounds, 0.0, true)
    nlp_raw = MathProgNLPModel(m)
    return run_IPOPT(nlp_raw, max_iter, bound_relax_factor,mehrotra_algorithm)
end

function run_IPOPT(nlp_raw, max_iter, bound_relax_factor, mehrotra_algorithm)
  t = 0
  hist = Array{ipopt_alg_history,1}()

  println(pd("t"), pd("comp"), pd("prm"), pd("dual"), pd("|x|"), pd("|y|"))
  A = jac(nlp_raw,nlp_raw.meta.x0)
  c = grad(nlp_raw,nlp_raw.meta.x0)

  ucon = nlp_raw.meta.ucon
  lcon = nlp_raw.meta.lcon
  #@show lbounds, ubounds
  infty = 1e30
  lbounds = max(-infty,nlp_raw.meta.lvar)
  ubounds = min(infty,nlp_raw.meta.uvar)
  #@show lbounds, ubounds

  for t = 1:max_iter

    solver = IpoptSolver(print_level=2, acceptable_iter=t, tol=1e-6, acceptable_tol=Inf, acceptable_compl_inf_tol=Inf, acceptable_constr_viol_tol=Inf, nlp_scaling_method="none", bound_relax_factor=bound_relax_factor, mehrotra_algorithm=mehrotra_algorithm)
    #setsolver(m, solver)
    #solve(m)
    #inner = m.internalModel.inner

    mp = NLPModels.NLPtoMPB(nlp_raw, solver)
    MathProgBase.optimize!(mp)
    inner = mp.inner


    hist_it = ipopt_alg_history()

    hist_it.t = t
    bounds_vio = max(0.0,maximum(lbounds - inner.x), max(0.0,maximum(inner.x -  ubounds )))
    a = cons(nlp_raw, inner.x)
    vio_A = max(maximum(a -  lcon), maximum(ucon - a))
    #norm(A * inner.x - b, Inf)

    hist_it.con_vio = max(vio_A, bounds_vio)
    #hist_it.primal_residual = max(norm(A * inner.x - b, Inf), max(0.0,maximum(lbounds - inner.x)), max(0.0,maximum(inner.x -  ubounds )))
    v1 = inner.mult_x_U - inner.mult_x_L
    v2 = A' * inner.mult_g
    v3 = v1 + v2


    hist_it.norm_grad_lag = norm(c + v3,Inf)
    hist_it.comp = max(maximum(abs(inner.x - lbounds) .* inner.mult_x_L),maximum(abs(inner.x - ubounds) .* inner.mult_x_U))
    hist_it.fval = dot(c, inner.x)
    hist_it.x_norm = norm(inner.x,Inf)
    hist_it.y_norm = max(norm(inner.mult_g,Inf),norm(inner.mult_x_U,Inf), norm(inner.mult_x_U,Inf))

    if t > 1 && hist_it.comp == hist[end].comp && hist_it.norm_grad_lag == hist[end].norm_grad_lag && hist_it.con_vio == hist[end].con_vio && hist_it.norm_grad_lag / (1.0 + norm(hist_it.y_norm,Inf)) < 1e-4 && hist[end].con_vio < 1e-4 && hist_it.comp < 1e-4
        break
    else
      println(pd(t), rd(hist_it.comp), rd(hist_it.con_vio), rd(hist_it.norm_grad_lag), rd(hist_it.x_norm), rd(hist_it.y_norm))
      push!(hist, hist_it)
    end
  end

  return hist
end


function load_netlib(num_nz::Int64)
    file_list = readdir("../netlib")
    problem_list = Array{String,1}()
    for file_name in file_list
        if file_name[(end-3):end] == ".mat" #&& filesize("netlib/$file_name") < file_size
          problem_name = file_name[1:(end-4)]
          println("loading...", problem_name)
          A, b, c, lbounds, ubounds = read_lp(problem_name)
          if nnz(A) <= num_nz
            push!(problem_list, problem_name)
          end
        end
    end

    return problem_list
end
