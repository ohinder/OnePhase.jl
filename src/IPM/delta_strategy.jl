function is_diag_dom(Q::SparseMatrixCSC{Float64,Int64})
    for i = 1:size(Q,2)
      if 3 * Q[i,i] < sum(Q[:,i]) + sum(Q[i,:])
        return false
      end
    end

    return true
end

function dist_diag_dom(Q::SparseMatrixCSC{Float64,Int64})
    delta_dist = 0.0
    for i = 1:size(Q,2)
      delta_dist = max(delta_dist,sum(Q[:,i]) + sum(Q[i,:]) - 3 * Q[i,i])
    end

    return delta_dist
end

#function bisection_search_delta!(iter::Class_iterate, kkt_solver::abstract_KKT_system_solver, min_delta::Float64, max_delta::Float64,  pars::Class_parameters, timer::class_advanced_timer)
#max_exponent = ceil(log(max_delta/min_delta))

function int_bisection(f::Function, a::Int64, b::Int64)
    for i = 1:10000
        if a == b - 1
          return b
        end
        midpoint = ceil((a + b)/2.0)
        if f(midpoint) == true
          b = midpoint
        else
          a = midpoint
        end
    end
end

function ipopt_strategy!(iter::Class_iterate, kkt_solver::abstract_KKT_system_solver, pars::Class_parameters, timer::class_advanced_timer)
    step_info = Blank_ls_info()

    MAX_IT = 500
    DELTA_ZERO = pars.delta.zero
    #DELTA_MIN = iter.point.mu * LinearAlgebra.norm(iter.point.y,Inf) / 100.0 #
    DELTA_MIN = pars.delta.min
    DELTA_MAX = pars.delta.max

    #Q = get_lag_hess(iter)
    #@show dist_diag_dom(Q)

    num_fac = 0

    inertia = 0
    status = :none

    tau = 1.5 * diag_min(kkt_solver)
    #Q = Symmetric(kkt_solver.M,:L)
    #eigvals,vec, = eigs(Q,nev=1,which=:SR,maxiter=10,tol=1e4)
    #@show eigvals
    #v = randn(length(iter.point.x))
    #v = v /LinearAlgebra.norm(v,2)
    #@show dot(v,Q * v)
    delta = DELTA_ZERO

    # see if we can succeed with delta = DELTA_ZERO
    if tau > 0.0
      tau = 0.0
      inertia = factor!( kkt_solver, delta, timer )
      num_fac += 1
      if inertia == 1
        return :success, num_fac, delta
      end
    end

    if pars.output_level >= 4
        println(pd("delta"), pd("inertia"))
    end

    for i = 1:MAX_IT
        if pars.output_level >= 4
          println("delta=",rd(delta), "inertia=", pd(inertia))
        end

        if i == 1
          if get_delta(iter) != 0.0
            delta = max(DELTA_MIN - tau, get_delta(iter) * pars.delta.dec)
          else
            delta = choose_delta_start(iter, kkt_solver, pars) - tau
          end
        else
            delta = delta * pars.delta.inc
        end

        inertia = factor!( kkt_solver, delta, timer )
        num_fac += 1

        n = length(iter.point.x)
        if inertia == 1
          return :success, num_fac, delta
        elseif is_diag_dom(kkt_solver.Q[1:n,1:n])
            println("WARNING: Inertia calculation incorrect")
            warn("Inertia calculation incorrect")
        end

        if delta > DELTA_MAX
          dx = LinearAlgebra.norm(kkt_solver.dir.x,Inf)
          dy = LinearAlgebra.norm(kkt_solver.dir.y,Inf)
          ds = LinearAlgebra.norm(kkt_solver.dir.s,Inf)

          my_warn("ipopt_strategy failed with delta_max=$DELTA_MAX, delta=$delta, i=$i")
          my_warn("num_fac=$num_fac, inertia=$inertia, status=$status, dir_x=$dx, dir_y=$dx, dir_s=$ds")

          #error("ipopt_strategy failed with too big a delta")
          return :failure, iter, delta
        end
    end

    error("max it")
end

function choose_delta_start(iter::Class_iterate, kkt_solver::abstract_KKT_system_solver, pars::Class_parameters)
    return pars.delta.start

    # do something clever to reduce number of factorizations

end
