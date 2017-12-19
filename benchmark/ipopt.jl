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
