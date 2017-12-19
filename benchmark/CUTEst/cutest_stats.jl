function filter_cutest(problem)
    correct_size = 50 <= problem["variables"]["number"] + problem["constraints"]["number"] && problem["constraints"]["number"] >= 10 && problem["variables"]["number"] <= 600 && problem["constraints"]["number"] <= 1000
    regular = problem["derivative_order"] >= 2 && problem["regular"] == true
    if correct_size && regular
        return true
    else
      return false
    end
end

#problem_list = CUTEst.select(max_var=100, max_var=1000, min_con=100, max_con=3000)
problem_list = CUTEst.select(custom_filter=filter_cutest)
problem_list = convert(Array{String,1},problem_list)

for problem_name in problem_list
  cute = CUTEstModel(problem_name)

  J = jac(cute, cute.meta.x0)
  H = hess(cute, cute.meta.x0, y=ones(size(J,2)))
  @show sum(H .> 0.0) / sum(J .> 0.0)
  finalize(cute)
end
