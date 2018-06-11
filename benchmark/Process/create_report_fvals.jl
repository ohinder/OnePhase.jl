include("../include.jl")
include("create_report.jl")
overlapping_results = get_CUTEst_results()

##
## FUNCTION VALUES
##
problem_list = collect(keys(first(overlapping_results)[2]))
method_list = collect(keys(overlapping_results))

CON_TOL = 1e-6
F_TOL = 0.5

winners = Dict()
winners["tie"] = 0
winners["fail"] = 0
begin
  for method_name in method_list
      print(pd(""), pd(method_name))
      winners[method_name] = 0
  end
  print("\n")

  for problem_name in problem_list
      print(pd(problem_name))
      best = "fail"
      both_kkt = true
      fbest = 1e14
      for method_name in method_list
        info = overlapping_results[method_name][problem_name]
        fval = info.fval
        con_vio = info.con_vio
        print(rd(con_vio))
        print(rd(fval))

        if info.status == :Optimal
          if con_vio > CON_TOL
            both_kkt = false
            println("CON VIOLATION!!")
          else
            ftol_scaled = F_TOL * (1.0 + abs(fbest))
            if fval <= fbest - ftol_scaled
              fbest = fval
              best = method_name
            elseif fbest - ftol_scaled < fval && fval <= fbest + ftol_scaled
              fbest = fval
              best = "tie"
            end
          end
        else
          both_kkt = false
        end
      end
      if both_kkt
        winners[best] += 1
      end
      print(pd(best))
      print("\n")
  end
end
