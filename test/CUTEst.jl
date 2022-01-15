using CUTEst
using CSV
optimal_problems = []
infeasible_problems = []
unbounded_problems = []
max_delta_problems = []
iteration_limit_problems = []
max_time_problems = []
error_problems = []

function cutest_tests()
    @testset "cutest" begin
        @testset "chain" begin
            nlp = CUTEstModel("CHAIN")
            try
				pars = OnePhase.Class_parameters()
		        pars.output_level = 0
	        	iter, status, hist, t, err, timer = OnePhase.one_phase_solve(nlp, pars);
                #@show OnePhase.get_original_x(iter) # gives the primal solution of the solver
                #@show OnePhase.get_y(iter) # gives the dual solution of the solver
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

function runModelFromProblem(cutest_problem, output_results_to_file)
    nlp = CUTEstModel(cutest_problem)
    try
        println("-----------EXECUTING PROBLEM----------", cutest_problem)
        pars = OnePhase.Class_parameters()
        pars.output_level = 0
        iter, status, hist, t, err, timer = OnePhase.one_phase_solve(nlp, pars);

		if output_results_to_file
			outputResultsToCSVFile(cutest_problem, hist)
		end

        @test true

		if status == :Optimal
			push!(optimal_problems, cutest_problem)
	    elseif status == :primal_infeasible
	        push!(infeasible_problems, cutest_problem)
	    elseif status == :dual_infeasible
	        push!(unbounded_problems, cutest_problem)
	    elseif status == :MAX_IT
	        push!(iteration_limit_problems, cutest_problem)
	    elseif status === :MAX_TIME
	        push!(max_time_problems, cutest_problem)
		else
			push!(max_delta_problems, cutest_problem)
	    end
    catch e
        push!(error_problems, cutest_problem)
        @test false
    finally
        finalize(nlp)
    end
end

function executeCUTEST_Models()
    problems = CUTEst.select(min_con = 1, max_con=5, only_free_var=true, max_var=5)

    for problem in problems
		try
	    	runModelFromProblem(problem, false)
    	catch e
	    	println(e)
    	end
    end

    println("------------RUNNING CUTEST TEST SET SUMMARY------------")
    println("----------------------------------------------Total: ", length(problems))
    println("--------------------------------------------OPTIMAL: ", length(optimal_problems))
    println("--------------------------------------------OPTIMAL: ", optimal_problems)
    println("-----------------------------------------INFEASIBLE: ", length(infeasible_problems))
    println("-----------------------------------------INFEASIBLE: ", infeasible_problems)
    println("------------------------------------------UNBOUNDED: ", length(unbounded_problems))
    println("------------------------------------------UNBOUNDED: ", unbounded_problems)
	println("------------------------------------ITERATION_LIMIT: ", length(iteration_limit_problems))
	println("------------------------------------ITERATION_LIMIT: ", iteration_limit_problems)
	println("-------------------------------------------MAX_TIME: ", length(max_time_problems))
	println("-------------------------------------------MAX_TIME: ", max_time_problems)
	println("------------------------------------------MAX_DELTA: ", length(max_delta_problems))
	println("------------------------------------------MAX_DELTA: ", max_delta_problems)
    println("----------------------------------------------ERROR: ", length(error_problems))
    println("----------------------------------------------ERROR: ", error_problems)
end


function executeCUTEST_Models_benchmark()
    #problems = CUTEst.select(min_con = 1, max_con=100, only_free_var=false, max_var=100)
    #problems = ["HIER13", "CHARDIS1", "HAGER2", "QPCBOEI1", "CHANDHEU", "MPC2", "BLOWEYA", "SINROSNB", "MPC10", "ACOPP57", "ACOPP300", "KISSING2", "LISWET9", "YAO", "VANDERM1", "LUBRIFC", "GMNCASE2", "CHEMRCTA", "CATENARY", "HYDROELM", "CATMIX", "ROTDISC", "LISWET7", "SSNLBEAM", "HAIFAL", "INTEGREQ", "MSQRTB", "DITTERT", "TRAINF", "ROBOTARM", "GPP", "YATP1SS", "QPNBOEI1", "LISWET6", "UBH5", "AGG", "HIER133A", "STEENBRG", "DRCAVTY1", "10FOLDTR", "EIGENA", "ARGLALE", "QR3DBD", "HIER16", "A0NSDSDS", "LEUVEN4", "OPTCDEG2", "YORKNET", "STEENBRA", "EIGMAXC", "MPC1", "TABLE5", "HADAMARD", "STEENBRC", "MPC11", "LEUVEN7", "BLOWEYB", "SOSQP2", "DTOC2", "DALLASM", "TABLE4", "TRO41X9", "CHAIN", "READING8", "QR3D", "QPCBOEI2", "KISSING", "SREADIN3", "LISWET4", "ACOPR118", "OPTCTRL6", "ORTHRGDS", "MPC8", "TRO21X5", "HANGING", "READING7", "CHANDHEQ", "FLOSP2HM", "DTOC4", "POLYGON", "MSS3", "FLOSP2TL", "GROUPING", "HIER133C", "LISWET10", "ORTHREGD", "EIGENACO", "ACOPR57", "FERRISDC", "DTOC1NC", "KSS", "TABLE3", "TABLE1", "VANDERM2", "CBRATU2D", "EIGENB2", "LEUVEN3", "EIGENB", "BROWNALE", "DTOC1ND", "DRCAVTY2", "LEUVEN5", "DTOC3", "SOSQP1", "DIXCHLNV", "STEENBRF", "CORKSCRW", "MSS2", "SPIN", "HAGER1", "MPC15", "EIGMINA", "C-RELOAD", "DRUGDISE", "ZAMB2", "KTMODEL", "HYDROELS", "ARWHDNE", "HAGER3", "LISWET3", "MSQRTA", "STEENBRE", "LEAKNET", "MOSARQP1", "LEUVEN6", "NGONE", "GMNCASE3", "MPC9", "FLOSP2HH", "DALLASL", "LINCONT", "FLOSP2HL", "TOYSARAH", "MPC5", "HIER163A", "CHEMRCTB", "JJTABEL3", "BLOWEYC", "BRATU3D", "ORTHRDS2", "BROYDNBD", "PRIMAL3", "TRAINH", "HVYCRASH", "EIGENBCO", "LEUVEN1", "GMNCASE4", "HIER163E", "DTOC1L", "SPINOP", "CBRATU3D", "EIGMINC", "LISWET12", "YATP1NE", "MOSARQP2", "DTOC1NA", "LISWET8", "MANNE", "ROSEPETAL", "A4X12", "GMNCASE1", "EIGMAXA", "HAGER4", "EIGMAXB", "WOODSNE", "HIE1327D", "HIE1372D", "ROCKET", "LISWET11", "ACOPR300", "HIER133E", "STEERING", "MADSSCHJ", "LISWET5", "ARGLCLE", "VANDERM3", "QPNBOEI2", "EIGENA2", "MPC3", "BROYDN3D", "OPTCDEG3", "LEUVEN2", "EIGENC", "SMMPSF", "EIGENC2", "STEENBRD", "NUFFIELD", "ACOPP118", "ARGTRIG", "LISWET2", "ELEC", "MPC13", "A5NSSNSM", "MPC12", "MINC44", "HIER163D", "LISWET1", "TABLE6", "HIER163C", "HIER133D", "READING3", "EIGENCCO", "DRCAVTY3", "A5NSDSDM", "LUBRIF", "CATENA", "DTOC1NB", "FLOSP2TM", "TABLE7", "MINPERM", "FLOSP2TH", "MPC14", "ZIGZAG", "SVANBERG", "SAWPATH", "ARGLBLE", "SPIN2", "QPNSTAIR", "EIGENAU", "QPCSTAIR", "HELSBY", "HIER133B", "MPC4", "POWELL20", "MPC7", "HIER163B", "OPTCTRL3", "HYDROELL", "ORTHREGF", "STEENBRB", "COSHFUN", "CAMSHAPE", "MPC16", "READING1", "EIGMINB", "SPIN2OP", "MPC6", "ORTHREGC"]
    problems = ["HIER13", "CHARDIS1", "HAGER2", "QPCBOEI1", "CHANDHEU", "MPC2", "BLOWEYA", "SINROSNB", "MPC10", "ACOPP57", "ACOPP300", "KISSING2", "LISWET9", "YAO", "VANDERM1", "LUBRIFC", "GMNCASE2", "CHEMRCTA", "CATENARY", "HYDROELM", "CATMIX", "ROTDISC", "LISWET7", "SSNLBEAM", "HAIFAL", "INTEGREQ", "MSQRTB", "DITTERT", "TRAINF", "ROBOTARM", "GPP", "QPNBOEI1", "LISWET6", "UBH5", "AGG", "HIER133A", "STEENBRG", "DRCAVTY1", "10FOLDTR", "EIGENA", "ARGLALE", "QR3DBD", "HIER16", "A0NSDSDS", "LEUVEN4", "OPTCDEG2", "YORKNET", "STEENBRA", "EIGMAXC", "MPC1", "TABLE5", "HADAMARD", "STEENBRC", "MPC11", "LEUVEN7", "BLOWEYB", "SOSQP2", "DTOC2", "DALLASM", "TABLE4", "TRO41X9", "CHAIN", "READING8", "QR3D", "QPCBOEI2", "KISSING", "SREADIN3", "LISWET4", "ACOPR118", "OPTCTRL6", "ORTHRGDS", "MPC8", "TRO21X5", "HANGING", "READING7", "CHANDHEQ", "FLOSP2HM", "DTOC4", "POLYGON", "MSS3", "FLOSP2TL", "GROUPING", "HIER133C", "LISWET10", "ORTHREGD", "EIGENACO", "ACOPR57", "FERRISDC", "DTOC1NC", "KSS", "TABLE3", "TABLE1", "VANDERM2", "CBRATU2D", "EIGENB2", "LEUVEN3", "EIGENB", "BROWNALE", "DTOC1ND", "DRCAVTY2", "LEUVEN5", "DTOC3", "SOSQP1", "DIXCHLNV", "STEENBRF", "CORKSCRW", "MSS2", "SPIN", "HAGER1", "MPC15", "EIGMINA", "C-RELOAD", "DRUGDISE", "ZAMB2", "KTMODEL", "HYDROELS", "ARWHDNE", "HAGER3", "LISWET3", "MSQRTA", "STEENBRE", "LEAKNET", "MOSARQP1", "LEUVEN6", "NGONE", "GMNCASE3", "MPC9", "FLOSP2HH", "DALLASL", "LINCONT", "FLOSP2HL", "TOYSARAH", "MPC5", "HIER163A", "CHEMRCTB", "JJTABEL3", "BLOWEYC", "BRATU3D", "ORTHRDS2", "BROYDNBD", "PRIMAL3", "TRAINH", "HVYCRASH", "EIGENBCO", "LEUVEN1", "GMNCASE4", "HIER163E", "DTOC1L", "SPINOP", "CBRATU3D", "EIGMINC", "LISWET12", "YATP1NE", "MOSARQP2", "DTOC1NA", "LISWET8", "MANNE", "ROSEPETAL", "A4X12", "GMNCASE1", "EIGMAXA", "HAGER4", "EIGMAXB", "WOODSNE", "HIE1327D", "HIE1372D", "ROCKET", "LISWET11", "ACOPR300", "HIER133E", "STEERING", "MADSSCHJ", "LISWET5", "ARGLCLE", "VANDERM3", "QPNBOEI2", "EIGENA2", "MPC3", "BROYDN3D", "OPTCDEG3", "LEUVEN2", "EIGENC", "SMMPSF", "EIGENC2", "STEENBRD", "NUFFIELD", "ACOPP118", "ARGTRIG", "LISWET2", "ELEC", "MPC13", "A5NSSNSM", "MPC12", "MINC44", "HIER163D", "LISWET1", "TABLE6", "HIER163C", "HIER133D", "READING3", "EIGENCCO", "DRCAVTY3", "A5NSDSDM", "LUBRIF", "CATENA", "DTOC1NB", "FLOSP2TM", "TABLE7", "MINPERM", "FLOSP2TH", "MPC14", "ZIGZAG", "SVANBERG", "SAWPATH", "ARGLBLE", "SPIN2", "QPNSTAIR", "EIGENAU", "QPCSTAIR", "HELSBY", "HIER133B", "MPC4", "POWELL20", "MPC7", "HIER163B", "OPTCTRL3", "HYDROELL", "ORTHREGF", "STEENBRB", "COSHFUN", "CAMSHAPE", "MPC16", "READING1", "EIGMINB", "SPIN2OP", "MPC6", "ORTHREGC"]
    for problem in problems
        runModelFromProblem(problem, true)
    end
	println("------------RUNNING CUTEST TEST SET SUMMARY------------")
    println("----------------------------------------------Total: ", length(problems))
    println("--------------------------------------------OPTIMAL: ", length(optimal_problems))
    println("--------------------------------------------OPTIMAL: ", optimal_problems)
    println("-----------------------------------------INFEASIBLE: ", length(infeasible_problems))
    println("-----------------------------------------INFEASIBLE: ", infeasible_problems)
    println("------------------------------------------UNBOUNDED: ", length(unbounded_problems))
    println("------------------------------------------UNBOUNDED: ", unbounded_problems)
	println("------------------------------------ITERATION_LIMIT: ", length(iteration_limit_problems))
	println("------------------------------------ITERATION_LIMIT: ", iteration_limit_problems)
	println("-------------------------------------------MAX_TIME: ", length(max_time_problems))
	println("-------------------------------------------MAX_TIME: ", max_time_problems)
	println("------------------------------------------MAX_DELTA: ", length(max_delta_problems))
	println("------------------------------------------MAX_DELTA: ", max_delta_problems)
    println("----------------------------------------------ERROR: ", length(error_problems))
    println("----------------------------------------------ERROR: ", error_problems)

    outputTotalResultsToFile(problems, optimal_problems, infeasible_problems, unbounded_problems, error_problems)
end

function outputResultsToCSVFile(cutest_problem, hist)
    current_dir = pwd()
    CSV.write("$current_dir/cutest_outputs/$cutest_problem.csv", hist, header = true)
end

function outputTotalResultsToFile(problems, optimal_problems, infeasible_problems, unbounded_problems, error_problems)
    current_dir = pwd()
    total = length(problems)
    total_optimal = length(optimal_problems)
    total_infeasible = length(infeasible_problems)
    total_unbounded = length(unbounded_problems)
    total_iteration_limit = length(iteration_limit_problems)
	total_max_time = length(max_time_problems)
	total_max_delta = length(max_delta_problems)
	total_erros = length(error_problems)

    open("$current_dir/cutest_outputs/total_results.txt" , "w") do file
	write(file, "----------------------------------------------Total: $total\n")
	write(file, "--------------------------------------------OPTIMAL: $total_optimal\n")
	write(file, "--------------------------------------------OPTIMAL: $optimal_problems\n")
	write(file, "-----------------------------------------INFEASIBLE: $total_infeasible\n")
	write(file, "-----------------------------------------INFEASIBLE: $infeasible_problems\n")
	write(file, "------------------------------------------UNBOUNDED: $total_unbounded\n")
	write(file, "------------------------------------------UNBOUNDED: $unbounded_problems\n")
	write(file, "------------------------------------ITERATION_LIMIT: $total_iteration_limitn\n")
	write(file, "------------------------------------ITERATION_LIMIT: $iteration_limit_problems\n")
	write(file, "-------------------------------------------MAX_TIME: $total_max_time\n")
	write(file, "-------------------------------------------MAX_TIME: $max_time_problems\n")
	write(file, "------------------------------------------MAX_DELTA: $total_max_delta")
	write(file, "------------------------------------------MAX_DELTA: $max_delta_problems\n")
	write(file, "----------------------------------------------ERROR: $total_erros\n")
	write(file, "----------------------------------------------ERROR: $error_problems")
    end
end
