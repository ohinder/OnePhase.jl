#=
0: Optimal solution found
-1: Time limit exceeded
-2: NAN occured
1: Iteration limit reached
2: Convergence to an infeasible point
4: Current solution estimate cannot be improved
5: Current solution estimate cannot be improved. It appears to be optimal,
but desired accuracy could not be achieved.

0: Optimal solution found
-1: Time limit exceeded
1: Maximal number of iterations exceeded
11: NaN or Inf occured
12: Number of iterations exceeded in restoration phase
16: Point is (almost) feasible, but restoration phase is called
17: Convergence to stationary point for infeasibility
18: Restoration phase cannot further improve feasibility
=#
function real_problems(problem)
    if problem["origin"] == "real" || problem["origin"] == "modelling"
        return true
    else
        return false
    end
end

function lrg_problem_func(min_num::Int64, max_num::Int64)
    function func(problem)
        regular = problem["derivative_order"] >= 2 && problem["regular"] == true && problem["constraints"]["number"] >= 1
        min_size_ok = problem["variables"]["number"] >= min_num && problem["constraints"]["number"] >= min_num
        #max_size = problem["variables"]["number"] <= max_num && problem["constraints"]["number"] <= max_num
        max_size_ok = problem["constraints"]["number"] + problem["variables"]["number"] <= max_num

        if min_size_ok && max_size_ok && regular
            return true
        else
          return false
        end
    end

    return func
end

function jump_status_conversion(info::problem_summary2; MAX_IT::Int64=3000)
  status = info.status
  if status == :Optimal || status == :first_order
    return :Optimal
  elseif status == :Infeasible || status == :infeasible
    return :primal_infeasible
  elseif status == :Unbounded
    return :dual_infeasible
  elseif status == :UserLimit
    if info.it_count >= MAX_IT
      return :MAX_IT
    else
      return :MAX_TIME
    end
  elseif status == :max_time
    return :MAX_TIME
  elseif status == :max_iter
    return :MAX_IT
  elseif status == :Error
    if info.it_count > 0
      return :ERROR
    else
      return :INIT_ERROR
    end
  else
    return :UNKNOWN
  end
end

function convert_JuMP(results::Dict{String,problem_summary2})
    for (prb_name, info) in results
      info.status = jump_status_conversion(info)
    end

    return results
end

function status_conversion(status::Int64)
    if status == 0
      return :Optimal
    elseif status == 17 || status == 2
      return :primal_infeasible
    #elseif status ==
    #  return :dual_infeasible
    elseif status == 1
      return :MAX_IT
    elseif status == -1
      return :MAX_TIME
    else
      return :ERROR
    end
end

function parse_number(something)
  if isa(something,Number)
      return something
  else
      return -1
  end
end

function load_results(file_name::String)
    results = readdlm(file_name,' ')
    results_dict = Dict{String, problem_summary2}()
    for i = 1:size(results,1)
      problem_name = results[i,1]
      results_dict[problem_name] = problem_summary2()

      results_dict[problem_name].it_count = parse_number(results[i,4])
      results_dict[problem_name].total_time = parse_number(results[i,7])
      results_dict[problem_name].status = status_conversion(parse_number(results[i,10]))
    end
    return results_dict
end

Infty = 99999999


function alg_success(status::Symbol)
    return (status == :Optimal  || status == :primal_infeasible || status == :dual_infeasible)
end


function iteration_list(results::Dict{String, problem_summary2}; MAX_IT::Int64=3000)
  problem_list = collect(keys(results))
  its = zeros(length(problem_list))

  for i = 1:length(problem_list)
      problem_name = problem_list[i]
      info = results[problem_name]

      if alg_success(info.status) && info.it_count < MAX_IT
        if info.it_count < 1
          println("warning: non-positive iteration count for $problem_name")
          its[i] = 1
        else
          its[i] = info.it_count
        end
      else
        its[i] = Infty
      end
  end

  return its
end

function best_its(it_list_by_method::Dict)
    n = length(first(values(it_list_by_method)))
    best_its = ones(n) * Infty
    for (method_name, it_list) in it_list_by_method
        @assert( length(it_list) == n )
        for i = 1:n
          best_its[i] = min(best_its[i], it_list[i])
        end
    end

    return best_its
end

function overlap(dics::Dict{String, Dict{String,problem_summary2}})
  overlap = collect(keys(first(dics)[2]))
  for (method_name, data) in dics
    problem_list = collect(keys(data))
    overlap = intersect(problem_list, overlap)
  end

  return restrict_problems(dics, overlap)

  #@show overlap
  #=new_dic = Dict{String, Dict{String,problem_summary2}}()
  for (method_name, data) in dics
    new_dic[method_name] = Dict{String,problem_summary2}()
    for problem_name in overlap
      new_dic[method_name][problem_name] = data[problem_name]
    end
  end=#
end


function all_status(res::Dict{String, Dict{String,problem_summary2}},status::Symbol)
  method_list = keys(res)

  opt_res = Dict{String, Dict{String,problem_summary2}}()
  for method in method_list
      opt_res[method] = Dict{String,problem_summary2}()
  end

  problem_list = keys(first(res)[2])
  for problem_name in problem_list
      include_prob = true
      for method in method_list
        if res[method][problem_name].status != status
          include_prob = false
        end
      end

      if include_prob
        for method in method_list
           opt_res[method][problem_name] = res[method][problem_name]
        end
      end
  end

  return opt_res
end

function restrict_to_set(dics::Dict{String, Dict{String,problem_summary2}}, set)
    p_list = []
    for (method_name,data) in dics
      p_list = union(get_list_of(data, set), p_list)
    end
    return restrict_problems(dics, p_list)
end


function get_list_of(results::Dict{String, problem_summary2}, set)
  new_problem_list = []
  for (problem_name, info) in results
      if info.status in set
        push!(new_problem_list,problem_name)
      end
  end

  return new_problem_list
end


function restrict_problems(results, problem_list)
  new_results = Dict{String,Dict{String,problem_summary2}}()
  for (method_name, sum_data) in results
      new_results[method_name] = Dict()
      for (problem_name, info) in sum_data
        if problem_name in problem_list
          new_results[method_name][problem_name] = info
        end
      end
  end
  return new_results
end

function tot(results::Dict{String, problem_summary2}, set)
  problem_list = collect(keys(results))
  tot = 0

  for i = 1:length(problem_list)
      problem_name = problem_list[i]
      info = results[problem_name]
      if info.status in set
        tot += 1
      end
  end

  return tot
end

function tot_failures(results::Dict{String, problem_summary2};MAX_IT::Int64=3000)
  problem_list = collect(keys(results))
  tot = 0

  for i = 1:length(problem_list)
      problem_name = problem_list[i]
      info = results[problem_name]
      if alg_success(info.status) && info.it_count <= MAX_IT
        # optimal
      else
        tot += 1
      end
  end

  return tot
end

function failure_list(results::Dict{String, problem_summary2};MAX_IT=3000)
  problem_list = collect(keys(results))
  failure_list = zeros(length(problem_list))

  for i = 1:length(problem_list)
      problem_name = problem_list[i]
      info = results[problem_name]
      if alg_success(info.status) && info.it_count <= MAX_IT
        failure_list[i] = 0
      else
        failure_list[i] = 1
      end
  end

  return failure_list
end

function shared_failures(results::Dict{String, Dict{String,problem_summary2}})
  problem_list = collect(keys(first(results)[2]))
  all_fail = ones(length(problem_list))
  #@show length(all_fail)
  for (method_name, data) in results
      ls = failure_list(data)
      #@show length(ls)
      all_fail = all_fail .* ls
  end
  return sum(all_fail)
end

function failure_causes(results::Dict{String, problem_summary2})
  problem_list = collect(keys(results))
  failure_causes = Dict()

  for i = 1:length(problem_list)
      problem_name = problem_list[i]
      info = results[problem_name]
      if !alg_success(info.status) #&& info.it_count <= 3000
        if info.status in keys(failure_causes)
          failure_causes[info.status] += 1
        else
          failure_causes[info.status] = 1
        end
      end
  end

  return failure_causes
end

function list_combined_failures(results::Dict{String, Dict{String,problem_summary2}})
  problem_list = collect(keys(first(results)[2]))
  all_fail = ones(length(problem_list))
  #@show length(all_fail)
  for (method_name, data) in results
      for (problem_name,info) in data
        ls = failure_list(data)
        #@show length(ls)
        all_fail = all_fail .* ls
      end
  end
  println("shared failures = ", sum(all_fail))
end

function list_combined_success(results::Dict{String, Dict{String,problem_summary2}})
  problem_list = collect(keys(first(results)[2]))
  all_succeed = ones(length(problem_list))
  #@show length(all_fail)
  for (method_name, data) in results
      for (problem_name,info) in data
        ls = 1.0 - failure_list(data)
        #@show length(ls)
        all_succeed = all_succeed .* ls
      end
  end
  println("shared successes = ", sum(all_succeed))
end

#=
function outcomes_table(results::Dict{String, Dict{String,problem_summary2}})
    outcomes = Dict()
    for problem_name in problem_list
        key = ()
        for (method_name, data) in results
            key = info.status
            #outcomes[(info.status,] += 1
        end
    end
end
=#



function print_failure_problems(results::Dict{String, Dict{String,problem_summary2}};MAX_IT::Int64=3000)
  problem_list = collect(keys(first(results)[2]))

  println("FAILURES")
  for i = 1:length(problem_list)
      problem_name = problem_list[i]
      first_failure = true

      for (method_name, data) in results
        info = results[method_name][problem_name]
        if alg_success(info.status) && info.it_count <= MAX_IT
          # optimal
        else
          if first_failure
            print(problem_name)
            first_failure = false
          end
          print(" ", method_name)
        end
      end
      if first_failure == false
        print("\n")
      end
  end
end

function remove_errors(dics::Dict{String, Dict{String,problem_summary2}}, error_set)
    new_dic = Dict{String, Dict{String,problem_summary2}}()
    for (method_name, data) in dics
      new_dic[method_name] = Dict{String,problem_summary2}()
      for (problem_name, info) in data
        if !(info.status in error_set)
          new_dic[method_name][problem_name] = info
        end
      end
    end
    return new_dic
end

function select(dics::Dict{String, Dict{String,problem_summary2}}, select_set)
    new_dic = Dict{String, Dict{String,problem_summary2}}()
    for (method_name, data) in dics
      new_dic[method_name] = Dict{String,problem_summary2}()
      for (problem_name, info) in data
        if (info.status in select_set)
          new_dic[method_name][problem_name] = info
        end
      end
    end
    return new_dic
end

function union_results(dics::Dict{String, Dict{String,problem_summary2}}, master::Dict{String, Dict{String,problem_summary2}})
    names_list = []
    for (method_name, data) in dics
      names_list = union(collect(keys(data)),names_list)
    end

    new_dic = Dict{String, Dict{String,problem_summary2}}()
    for (method_name, data) in results
      new_dic[method_name] = Dict{String,problem_summary2}()
      for (problem_name, info) in data
        if problem_name in names_list
          new_dic[method_name][problem_name] = info
        end
      end
    end
    return new_dic
end

function show_results(results::Dict{String, Dict{String,problem_summary2}})
    problem_list = collect(keys(first(results)[2]))
    for problem_name in problem_list
      for (method_name, data) in results
         info = data[problem_name]
         println(method_name, " ", problem_name, " ", info.it_count, " ", info.status)
      end
    end
end


function compute_its_etc(overlapping_results; MAX_IT::Int64=3000)
    its = Dict{String,Array{Int64,1}}()
    println("Quartiles of iterations:")
    for (method_name, sum_data) in overlapping_results
        its[method_name] = iteration_list(sum_data, MAX_IT=MAX_IT);
        d = its[method_name]
        println(OnePhase.pd(method_name,20), " = ", OnePhase.rd(quantile(d,0.25)), OnePhase.rd(quantile(d,0.5)), OnePhase.rd(quantile(d,0.75)))
    end

    println("Times:")
    times = Dict{String,Array{Float64,1}}()
    for (method_name, sum_data) in overlapping_results
        times[method_name] = []
        for (problem_name,info) in sum_data
          if alg_success(info.status)
              t = info.total_time
          else
              t = Infty
          end

          push!(times[method_name], t);
        end
        tm = times[method_name]
        println(OnePhase.pd(method_name,20), " = ", OnePhase.rd(quantile(tm,0.25)), OnePhase.rd(quantile(tm,0.5)), OnePhase.rd(quantile(tm,0.75)))
        #println(method_name, " mean=", rd(mean(times[method_name])), " median=", rd(median(times[method_name])))
    end

    println("Quartiles of iteration ratios to best solver:")
    best = best_its(its)

    ratios = Dict()
    for (method_name, val) in its
      ratios[method_name] = its[method_name] ./ best;
      lrg = its[method_name] .>= Infty
      ratios[method_name][lrg] = Inf

      d = ratios[method_name]
      println(OnePhase.pd(method_name,20), " = ", OnePhase.rd(quantile(d,0.25)), OnePhase.rd(quantile(d,0.5)), OnePhase.rd(quantile(d,0.75)))
    end
    return its, best, ratios, times
end
####################################################################################################
####################################################################################################
#############################       DUAL VARIABLE STUFF                #############################
####################################################################################################
####################################################################################################


function compute_dual_hist(results)
    dual_hist = Array{Float64,1}()
    t = 0
    for it_hist in results
      if it_hist.t > t
        push!( dual_hist, it_hist.y_norm )
        t += 1
      end
    end
    #dual_maxes[folder_name][problem_name] = maximum(dual_hist)

    return dual_hist
end

function get_lists(summary)
    problem_list = collect(keys(first(summary)[2]))
    method_list = keys(summary)

    return problem_list, method_list
end

function get_all_optimal(summary, problem_list) # get all problems where all solvers get the optimal solution
    remove_these = Array{String,1}()

    for problem_name in problem_list
      for (method_name, sum_data) in summary_netlib
          if sum_data[problem_name].status != :Optimal
            push!(remove_these, problem_name)
          end
      end
    end

    all_optimal = Array{String,1}()

    for problem_name in problem_list
        if !(problem_name in remove_these)
          push!(all_optimal, problem_name)
        end
    end

    return all_optimal
end

function summarize_by_iteration(hist)
    sum_hist = Array{abstract_alg_history,1}()
    for i = 2:length(hist)
      if hist[i].t > hist[i-1].t
        push!(sum_hist, hist[i-1])
      end
    end
    push!(sum_hist, hist[end])

    return sum_hist
end


function aggregate_summary_by_it_for_plots(sh, field_list::Dict{String,Symbol})
    data_us = Dict()

    for (label, field) in field_list
        data_us[field] = zeros(length(sh));
    end

    for i = 1:length(sh)
        for (label, field) in field_list
          data_us[field][i] = getfield(sh[i],field);
        end
    end

    return data_us
end
