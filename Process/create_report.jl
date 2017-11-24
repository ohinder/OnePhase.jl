function get_CUTEst_results()
  results = Dict{String, Dict{String,problem_summary}}()

  results["one phase"] = load("../results/one_phase/sept_3_corrections/summary.jld", "summary")
  results["ipopt"] = convert_JuMP(load("../results/ipopt/plain/summary.jld", "summary"))

  min_num = 100
  max_num = 10000
  lrg_problems = lrg_problem_func(min_num, max_num)
  problem_list = CUTEst.select(custom_filter=lrg_problems);
  problem_list = convert(Array{String,1},problem_list);

  results = restrict_problems(results, problem_list)

  return results
end
