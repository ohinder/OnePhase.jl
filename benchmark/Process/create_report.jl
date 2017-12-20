
#folder = "CUTEst_infeasible"
#folder = "CUTEst"
folder = "CUTEst_low_tol"
#folder = "netlib_infeasible"

if folder == "CUTEst"
  mode = :optimal
  data = :CUTEst
elseif folder == "CUTEst_low_tol"
  mode = :optimal
  data = :CUTEst_low_tol
elseif folder == "CUTEst_infeasible"
  mode = :primal_infeasible
  data = :CUTEst
elseif folder == "netlib_infeasible"
  mode = :primal_infeasible
  data = :netlib
end


function get_CUTEst_results()
  results = Dict{String, Dict{String,problem_summary2}}()

  #GLOBAL mode;
  if mode == :optimal
    if data == :CUTEst
    #results["ipopt"] = load("../results/one_phase/sept_3_corrections/summary.jld", "summary")
      results["one phase"] = cps(load("../results/one_phase/Dec1/summary.jld", "summary"))
      results["ipopt"] = convert_JuMP(cps(load("../results/ipopt/plain/summary.jld", "summary")))
    elseif data == :CUTEst_low_tol
      results["one phase"] = cps(load("../results/one_phase/LowTol_Dec19/summary.jld", "summary"))
      results["ipopt"] = convert_JuMP(cps(load("../results/ipopt/very_low_tol/summary.jld", "summary")))
    end
  elseif mode == :primal_infeasible
    if data == :CUTEst
      results["one phase"] = cps(load("../results/one_phase/infeas-Dec11/summary.jld", "summary"))
      results["ipopt"] = convert_JuMP(cps(load("../results/ipopt/infeas-Dec15/summary.jld", "summary")))
    elseif data == :netlib
      results["one phase"] = cps(load("../results/one_phase/infeas-lp-test/summary.jld", "summary"))
      results["ipopt"] = convert_JuMP(cps(load("../results/ipopt/infeas-lp-test/summary.jld", "summary")))
    end
  else
    error("")
  end

  #=min_num = 0
  max_num = 10000
  lrg_problems = lrg_problem_func(min_num, max_num)
  problem_list = CUTEst.select(custom_filter=lrg_problems);
  problem_list = convert(Array{String,1},problem_list);

  results = restrict_problems(results, problem_list)=#

  overlapping_results = overlap(results)

  return overlapping_results
end
