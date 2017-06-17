include("include.jl")

using MAT

function read_lp(name::String)
    lp_data = matopen("netlib/$(name).mat")

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

using Ipopt, JuMP
function read_lp_into_JuMP(name::String, solver)
  A, b, c, lbounds, ubounds = read_lp(name)

  threshold = 1e20
  li = (lbounds .> -threshold)
  ui = (ubounds .< threshold)

  mod = Model(solver=solver)

  n = length(c)
  x = Array{JuMP.Variable,1}(n)
  for i = 1:n
    if true == li[i] && true == ui[i]
      x[i] = @variable(mod, lowerbound = lbounds[i], upperbound = ubounds[i])
    elseif true == li[i]
      x[i] = @variable(mod, lowerbound = lbounds[i])
    elseif true == ui[i]
      x[i] = @variable(mod, upperbound = ubounds[i])
    else
      x[i] = @variable(mod)
    end
  end

  @objective(mod, Min, (c' * x)[1])
  @constraint(mod, A * x .== b)

  return mod
end


function run_netlib_problems_using_IPOPT(problems::Array{String,1}, test_name::String)
    summary = Dict{String, problem_summary}()

    folder_name = "results/$test_name"
    jld_folder = "$folder_name/jld"
    log_folder = "$folder_name/log"
    if_mkdir(folder_name)
    if_mkdir(jld_folder)
    if_mkdir(log_folder)

    for problem_name in problem_list
        start_time = time()
        println("RUNNING ", problem_name)

        ORG_STDOUT = STDOUT
        file = open("$log_folder/$(problem_name).txt", "w")
        redirect_stdout(file)

        A, b, c, lbounds, ubounds = read_lp(problem_name)

        hist = Array{ipopt_alg_history,1}()

        println(pd("t"), pd("comp"), pd("prm"), pd("dual"), pd("dx"), pd("dy"))
        status = :Blank
        t = 0
        for t = 1:100
          solver = IpoptSolver(print_level=2, acceptable_iter=t, acceptable_tol=Inf, acceptable_compl_inf_tol=Inf, acceptable_constr_viol_tol=Inf)
          m = read_lp_into_JuMP(problem_name, solver)
          status = solve(m)

          inner = m.internalModel.nlpmodel.inner

          hist_it = ipopt_alg_history()

          hist_it.t = t
          hist_it.primal_residual = max(norm(A * inner.x - b, Inf), max(0.0,maximum(lbounds - inner.x)), max(0.0,maximum(inner.x -  ubounds )))
          hist_it.norm_grad_lag = norm(c + A' * inner.mult_g + inner.mult_x_U - inner.mult_x_L,Inf)
          hist_it.comp = max(maximum(abs(inner.x - lbounds) .* inner.mult_x_L),maximum(abs(inner.x - ubounds) .* inner.mult_x_U))
          hist_it.fval = dot(c, inner.x)
          hist_it.x_norm = norm(inner.x,Inf)
          hist_it.y_norm = max(norm(inner.mult_g,Inf),norm(inner.mult_x_U,Inf), norm(inner.mult_x_U,Inf))



          if t > 1 && hist_it.norm_grad_lag == hist[end].norm_grad_lag && hist_it.primal_residual == hist[end].primal_residual && hist_it.norm_grad_lag < 1e-4 && hist[end].primal_residual < 1e-4
              break
          else
            println(pd(t), rd(hist_it.comp), rd(hist_it.primal_residual), rd(hist_it.norm_grad_lag), rd(hist_it.x_norm), rd(hist_it.y_norm))
            push!(hist, hist_it)
          end
        end

        solver = IpoptSolver(print_level=5, max_iter=100)
        m = read_lp_into_JuMP(problem_name, solver)
        status = solve(m)

        summary[problem_name] = problem_summary()
        summary[problem_name].status = status;
        summary[problem_name].total_time = time() - start_time
        set_info_me!(summary[problem_name], hist)

        redirect_stdout(STDOUT)
        close(file)

        save("$(folder_name)/summary.jld","summary",summary)
        save("$(jld_folder)/$(problem_name).jld","history",hist)

        summary_file = open("$(folder_name)/summary.txt", "w")
        write_summary(summary_file, summary)
        close(summary_file)
    end
end

function run_netlib_problems_using_our_solver(problems::Array{String,1}, test_name::String, par::Class_parameters)
    summary = Dict{String, problem_summary}()

    folder_name = "results/$test_name"
    jld_folder = "$folder_name/jld"
    log_folder = "$folder_name/log"
    if_mkdir(folder_name)
    if_mkdir(jld_folder)
    if_mkdir(log_folder)

    for problem_name in problem_list
        start_time = time()
        println("RUNNING ", problem_name)
        nlp = read_lp_into_Class_QP(problem_name);

        ORG_STDOUT = STDOUT
        file = open("$log_folder/$(problem_name).txt", "w")
        redirect_stdout(file)

        reset_advanced_timer()
        start_advanced_timer()
        intial_it = init(nlp, my_par);
        @assert(is_feasible(intial_it, my_par.comp_feas))
        iter, status, hist, t, err = one_phase_IPM(intial_it, my_par);
        pause_advanced_timer()
        print_timer_stats()

        save("$(jld_folder)/$(problem_name).jld","history",hist)
        summary[problem_name] = problem_summary()
        summary[problem_name].status = status
        set_info_me!(summary[problem_name], hist)
        summary[problem_name].total_time = time() - start_time

        redirect_stdout(ORG_STDOUT)
        close(file)

        save("$(folder_name)/summary.jld", "summary", summary, "pars", par)

        summary_file = open("$(folder_name)/summary.txt", "w")
        write_summary(summary_file, summary)
        close(summary_file)
    end
end

#nlp = read_lp("AGG");
#nlp = read_lp("AGG");
#problem_list = ["BANDM", "SC50A", "SC50B", "SC105", "ADLITTLE", "BORE3D", "STOCFOR1", "KB2", "SHARE2B", "AFIRO", "BRANDY"]#, "BNL1"] #,"BNL1","BORE3D","AGG"]
#problem_list = ["AGG"]
file_list = readdir("netlib")
problem_list = Array{String,1}()
for file_name in file_list
    if file_name[(end-3):end] == ".mat" && filesize("netlib/$file_name") < 5e4
      push!(problem_list,file_name[1:(end-4)])
    end
end

problem_list



if true
test_name = "netlib"
run_netlib_problems_using_our_solver(problem_list, test_name, my_par)
end

if false
test_name = "netlib-ipopt"
run_netlib_problems_using_IPOPT(problem_list, test_name)
end

if false
solver = IpoptSolver(print_level=7)
m = read_lp_into_JuMP("AGG3", solver)
solve(m)
end
