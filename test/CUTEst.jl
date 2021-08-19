using CUTEst

function cutest_tests()
    @testset "cutest" begin
        @testset "chain" begin
            nlp = CUTEstModel("CHAIN")
            try
                iter, status, hist, t, err, timer = OnePhase.one_phase_solve(nlp);
                @show OnePhase.get_original_x(iter) # gives the primal solution of the solver
                @show OnePhase.get_y(iter) # gives the dual solution of the solver
                @test true
            catch(e)
                println(e)
                @test false
            finally
                finalize(nlp)
            end
        end
    end
end
