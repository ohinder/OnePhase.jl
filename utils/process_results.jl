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

function jump_status_conversion(status::Symbol)
  if status == :Optimal
    return :optimal
  elseif status == :Infeasible
    return :primal_infeasible
  elseif status == :Unbounded
    return :dual_infeasible
  #elseif status ==
  #  return :dual_infeasible
elseif status == :UserLimit
    return :MAX_IT
  elseif status == -1
    return :MAX_TIME
  else
    return :ERROR
  end
end

function convert_JuMP(results::Dict{String,problem_summary})
    for (prb_name, info) in results
      info.status = jump_status_conversion(info.status)
    end

    return results
end

function status_conversion(status::Int64)
    if status == 0
      return :optimal
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
    results_dict = Dict{String, problem_summary}()
    for i = 1:size(results,1)
      problem_name = results[i,1]
      results_dict[problem_name] = problem_summary()

      results_dict[problem_name].it_count = parse_number(results[i,4])
      results_dict[problem_name].total_time = parse_number(results[i,7])
      results_dict[problem_name].status = status_conversion(parse_number(results[i,10]))
    end
    return results_dict
end

Infty = 99999999


function alg_success(status::Symbol)
    return (status == :optimal  || status == :primal_infeasible || status == :dual_infeasible)
end


function iteration_list(results::Dict{String, problem_summary})
  problem_list = collect(keys(results))
  its = zeros(length(problem_list))

  for i = 1:length(problem_list)
      problem_name = problem_list[i]
      info = results[problem_name]
      if alg_success(info.status) && info.it_count < 3000
        its[i] = info.it_count
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

function overlap(dics::Dict{String, Dict{String,problem_summary}})
  overlap = collect(keys(first(dics)[2]))
  for (method_name, data) in dics
    problem_list = collect(keys(data))
    overlap = intersect(problem_list, overlap)
  end

  #@show overlap
  new_dic = Dict{String, Dict{String,problem_summary}}()
  for (method_name, data) in dics
    new_dic[method_name] = Dict{String,problem_summary}()
    for problem_name in overlap
      new_dic[method_name][problem_name] = data[problem_name]
    end
  end

  return new_dic
end

function tot_failures(results::Dict{String, problem_summary})
  problem_list = collect(keys(results))
  tot = 0

  for i = 1:length(problem_list)
      problem_name = problem_list[i]
      info = results[problem_name]
      if alg_success(info.status) && info.it_count <= 3000
        # optimal
      else
        tot += 1
      end
  end

  return tot
end

function list_failures(results::Dict{String, Dict{String,problem_summary}})
  problem_list = collect(keys(first(results)[2]))

  println("FAILURES")
  for i = 1:length(problem_list)
      problem_name = problem_list[i]
      first_failure = true

      for (method_name, data) in results
        info = results[method_name][problem_name]
        if alg_success(info.status) && info.it_count <= 3000
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

function remove_errors(dics::Dict{String, Dict{String,problem_summary}}, error_set)
    new_dic = Dict{String, Dict{String,problem_summary}}()
    for (method_name, data) in dics
      new_dic[method_name] = Dict{String,problem_summary}()
      for (problem_name, info) in data
        if !(info.status in error_set)
          new_dic[method_name][problem_name] = info
        end
      end
    end
    return new_dic
end

function select(dics::Dict{String, Dict{String,problem_summary}}, select_set)
    new_dic = Dict{String, Dict{String,problem_summary}}()
    for (method_name, data) in dics
      new_dic[method_name] = Dict{String,problem_summary}()
      for (problem_name, info) in data
        if (info.status in select_set)
          new_dic[method_name][problem_name] = info
        end
      end
    end
    return new_dic
end

function union_results(dics::Dict{String, Dict{String,problem_summary}}, master::Dict{String, Dict{String,problem_summary}})
    names_list = []
    for (method_name, data) in dics
      names_list = union(collect(keys(data)),names_list)
    end

    new_dic = Dict{String, Dict{String,problem_summary}}()
    for (method_name, data) in results
      new_dic[method_name] = Dict{String,problem_summary}()
      for (problem_name, info) in data
        if problem_name in names_list
          new_dic[method_name][problem_name] = info
        end
      end
    end
    return new_dic
end

function show_results(results::Dict{String, Dict{String,problem_summary}})
    problem_list = collect(keys(first(results)[2]))
    for problem_name in problem_list
      for (method_name, data) in results
         info = data[problem_name]
         println(method_name, " ", problem_name, " ", info.it_count, " ", info.status)
      end
    end
end


function compute_its_etc(overlapping_results)
    its = Dict{String,Array{Int64,1}}()
    println("quartiles")
    for (method_name, sum_data) in overlapping_results
        its[method_name] = iteration_list(sum_data);
        d = its[method_name]
        println(pd(method_name,20), " = ", rd(quantile(d,0.2)), rd(quantile(d,0.4)), rd(quantile(d,0.6)), rd(quantile(d,0.8)))
    end

    times = Dict{String,Array{Float64,1}}()
    for (method_name, sum_data) in overlapping_results
        times[method_name] = []
        for (problem_name,info) in sum_data
          push!(times[method_name], info.total_time);
        end
        println(method_name, " ", rd(mean(times[method_name])), " ", rd(median(times[method_name])))
    end

    best = best_its(its)

    ratios = Dict()
    for (method_name, val) in its
      ratios[method_name] = its[method_name] ./ best;
      lrg = its[method_name] .>= Infty
      ratios[method_name][lrg] = Inf

      d = ratios[method_name]
      println(pd(method_name,20), " = ", rd(quantile(d,0.2)), rd(quantile(d,0.4)), rd(quantile(d,0.6)), rd(quantile(d,0.8)), rd(quantile(d,0.9)))
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
          if sum_data[problem_name].status != :optimal
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
