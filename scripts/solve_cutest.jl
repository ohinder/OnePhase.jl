import ArgParse
include("../benchmark/CUTest/run_cutest.jl")

"""
Defines parses and args.
# Returns
A dictionary with the values of the command-line arguments.
"""
function parse_command_line()
  arg_parse = ArgParse.ArgParseSettings()

  help_method = "The optimization method to use, must be `one-phase` or `ipopt`."

  ArgParse.@add_arg_table! arg_parse begin
    "--method"
    help = help_method
    arg_type = String
    required = true

    "--output_dir"
    help = "The directory for output files."
    arg_type = String
    required = true

    "--linear_solver"
    help = "The linear solver type to use, must be ':julia' or 'HSL'"
    arg_type = Symbol
    required = true

    "--kkt_solver"
    help = "The kkt solver type to use, when selecting one-phase it must be ':schur' or ':symmetric' or ':clever_symmetric'"
    arg_type = Symbol

    "--max_it"
    help = "The maximum number of iterations to run"
    arg_type = Int64
    default = 3000

    "--max_time"
    help = "The maximum time to run in seconds"
    arg_type = Float64
    default =  60.0^2

    "--tol_opt"
    help = "The tolerance for optimality"
    arg_type = Float64
    default =  1e-6

    "--verbosity"
    help =
      "The verbosity level for printing. Values between 1 and 5 " *
      "print generic information on the solve process, i.e., a table of the " *
      "iteration statistics. Values greater than 5 provide information " *
      "useful for developers."
    arg_type = Int64
    default = 2

    "--min_nvar"
    help = "The minimum number of variables for CUTEst model"
    arg_type = Int64
    default =  10

    "--max_nvar"
    help = "The maximum number of variables for CUTEst model"
    arg_type = Int64
    default =  10000

    "--min_ncon"
    help = "The minimum number of constraints for CUTEst model"
    arg_type = Int64
    default = 10

    "--max_ncon"
    help = "The maximum constraints of variables for CUTEst model"
    arg_type = Int64
    default =  10000

  end

  return ArgParse.parse_args(arg_parse)
end

function main()
  parsed_args = parse_command_line()

  if parsed_args["method"] == "one-phase"
    my_par = OnePhase.Class_parameters()
    my_par.term.max_time = 60.0 * 60
    my_par.term.max_it = 3000
    if parsed_args["linear_solver"] ==  Symbol(":julia") || parsed_args["linear_solver"] ==  Symbol(":HSL")
      my_par.kkt.linear_solver_type = parsed_args["linear_solver"]
      if parsed_args["linear_solver"] ==  Symbol(":HSL")
        OnePhase.setUSE_HSL(true)
        OnePhase.loadHSL("../src/linear_system_solvers/")
      end
    else
     error("Unknown linear solver type")
    end

    if parsed_args["kkt_solver"] ==  Symbol(":schur") || parsed_args["kkt_solver"] ==  Symbol(":symmetric") || parsed_args["kkt_solver"] ==  Symbol(":clever_symmetric")
      my_par.kkt.kkt_solver_type = parsed_args["kkt_solver"]
    else
     error("Unknown kkt solver type")
    end

    my_par.term.max_it = parsed_args["max_it"]
    my_par.term.max_time = parsed_args["max_time"]
    my_par.term.tol_opt = parsed_args["tol_opt"]
    my_par.output_level = parsed_args["verbosity"]
    min_ncon = parsed_args["min_ncon"]
    max_ncon = parsed_args["max_ncon"]
    min_nvar = parsed_args["min_nvar"]
    max_nvar = parsed_args["max_nvar"]

    folder_name = parsed_args["output_dir"]
    println("pwd: ", pwd())
    # if_mkdir("../benchmark/results/$folder_name")
    if_mkdir("$folder_name")
    run_cutest_problems_using_our_solver(folder_name, my_par, min_ncon, max_ncon, min_nvar, max_nvar)
  elseif parsed_args["method"] == "ipopt"
    linear_solver = ""
    if parsed_args["linear_solver"] ==  Symbol(":HSL")
      linear_solver = "ma97"
    end
    max_it = parsed_args["max_it"]
    max_time = parsed_args["max_time"]
    tol_opt = parsed_args["tol_opt"]
    output_level = parsed_args["verbosity"]
    min_ncon = parsed_args["min_ncon"]
    max_ncon = parsed_args["max_ncon"]
    min_nvar = parsed_args["min_nvar"]
    max_nvar = parsed_args["max_nvar"]

    folder_name = parsed_args["output_dir"]
    # if_mkdir("../benchmark/results/$folder_name")
    if_mkdir("$folder_name")
    run_cutest_problems_on_solver(folder_name, output_level, tol_opt, max_it, max_time, linear_solver, min_ncon, max_ncon, min_nvar, max_nvar)
  else
    error("`method` arg must be either `one-phase` or `ipopt`.")
  end
end

main()
