include("run_cutest.jl")
using Ipopt
TIME_LIMIT = 60.0^2
problem_list = default_list()

if false
  folder_name = "ipopt/very_low_tol"
  if_mkdir("../results/$folder_name")
  run_cutest_problems_on_solver(problem_list, folder_name, print_level=3, tol=1e-2, max_iter=3000, max_cpu_time=TIME_LIMIT, nlp_scaling_method="none", bound_relax_factor = 0.0, acceptable_iter=999999)
end

if true
  folder_name = "ipopt/Jan2"
  if_mkdir("../results/$folder_name")
  run_cutest_problems_on_solver(problem_list, folder_name, 4, 1e-6, 3000, TIME_LIMIT, "none", 0.0, 999999)
end

if false
  folder_name = "ipopt/medium_perturb"
  if_mkdir("../results/$folder_name")
  run_cutest_problems_on_solver(problem_list, folder_name, print_level=3, tol=1e-6, max_iter=3000, max_cpu_time=TIME_LIMIT, bound_relax_factor = 1e-8, nlp_scaling_method="none", acceptable_iter=Inf)
end

if false
  folder_name = "ipopt/medium_high_tol"
  if_mkdir("../results/$folder_name")
  run_cutest_problems_on_solver(problem_list, folder_name, print_level=3, tol=1e-8, max_iter=3000, max_cpu_time=TIME_LIMIT, bound_relax_factor = 0.0, nlp_scaling_method="none", acceptable_iter=Inf)
end

if false
  folder_name = "ipopt/medium_perturb_high_tol"
  if_mkdir("../results/$folder_name")
  run_cutest_problems_on_solver(problem_list, folder_name, print_level=3, tol=1e-8, max_iter=3000, max_cpu_time=TIME_LIMIT, nlp_scaling_method="none", acceptable_iter=Inf)
end
