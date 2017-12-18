include("Class_point.jl")
include("Class_iterate.jl")
include("misc.jl")
include("eval.jl")
include("Class_cutest.jl")
include("summary.jl")
include("process_results.jl")
include("plots.jl")
include("autotune.jl")


function IpoptSolve(nlp)
    tic()
    solver = IpoptSolver(print_level=5, max_iter=3000, bound_relax_factor=0.0, nlp_scaling_method="none", acceptable_iter=999999)
    #,mehrotra_algorithm="yes") #, tol_dual_abs=1e-6)
    #solver = IpoptSolver(print_level=5, tol=1e-8)
    mp = NLPModels.NLPtoMPB(nlp, solver)
    MathProgBase.optimize!(mp)
    @show norm(mp.inner.mult_g, Inf)
    #y = MathProgBase.getdual(mp)
    solver = MathProgBase.getrawsolver(mp)
    #finalize(nlp_raw2)
    toc();
end


function get_file_names(dir::String,ext::String)
  name_list = readdir(dir)
  cleaned_name_list = []
  for name_full in name_list
    try
      if name_full[(end-3):end] == ".mat"
          name = name_full[1:(end-4)]
      end
      cleaned_name_list = [cleaned_name_list; name]
      catch(e)
        @show e
        println("ERROR in running " * name_full)
      end
  end

  return cleaned_name_list
end
