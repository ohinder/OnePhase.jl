type woodbury_identity
    U::AbstractMatrix
    V::AbstractMatrix
    invA_U::AbstractMatrix
    factored_capacitance_matrix
    factored_A::abstract_linear_system_solver

    # (A + U*V')^-1
    function woodbury_identity(factored_A::abstract_linear_system_solver, U::AbstractArray, V::AbstractArray)
        invA_U = false
        capacitance_matrix = false
        try
          invA_U = ls_solve(factored_A, U)
          n = size(U,1);
          k = size(U,2);

          capacitance_matrix = (eye(k) + V * invA_U);
          factored_capacitance_matrix = lufact(capacitance_matrix);

          return new(U, V, invA_U, factored_capacitance_matrix, factored_A);
        catch e
            println("ERROR woodbury_identity")
            #@show full(factored_A._SparseMatrix)
            #@show capacitance_matrix
            #@show invA_U'
            throw(e)
        end
    end
end

function evaluate(W::woodbury_identity, sol::Array{Float64,1})
    return W.factored_A._SparseMatrix * sol + W.U * (W.V * sol)
end

function rel_error(W::woodbury_identity, rhs::Array{Float64,1}, sol::Array{Float64,1})
    return norm(evaluate(W, sol) - rhs,1)/norm(rhs,1);
end

function ls_solve(W::woodbury_identity, rhs::Array{Float64,1})
  try
    rhs_err = rhs
    sol = zeros(length(rhs))
    for i = 1:1
        # iterative refinement
        my_temp = ls_solve(W.factored_A, rhs_err)
        Δd = (my_temp - W.invA_U * (W.factored_capacitance_matrix \ (W.V * my_temp)))
        sol = sol + Δd[:]
        rhs_err = rhs - evaluate(W,sol)
        #@show i, norm(rhs_err,1)
    end

    # catch errors
    tol = 1e-8;
    err = rel_error(W, rhs, sol) # why is this rel error?
    if err > tol
        warn("numerical instability using Woodbury, using direct factorization instead of woodbury.")
        sol = ls_solve_direct(W, rhs)
        err = rel_error(W, rhs, sol)

        if (err > tol)
            error("numerical stability using direct factorization, computation skipped")
        end
    end

    return sol
  catch e
      println("ERROR ls_solve")
      throw(e)
  end
end

function ls_solve_direct(W::woodbury_identity, rhs::Array{Float64,1}) # if there is numerical errors
    mat = W.factored_A._SparseMatrix + sparse(W.U) * sparse(W.V)
    return lufact(mat) \ rhs;
end

function numerical_error(A, sol, rhs)
    return norm(A * sol - rhs,1);
end

function e_(i::Int64,n::Int64)
    @assert(i <= n)
    vec = zeros(n);
    vec[i] = 1.0;
    return vec;
end
