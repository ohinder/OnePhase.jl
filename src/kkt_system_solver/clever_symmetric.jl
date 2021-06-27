using LinearAlgebra
# UNDER CONSTRUCTION

mutable struct Parallel_row
    ind::Int64
    ratio::Float64
    u::Float64
    g::Float64

    function Parallel_row(ind::Int64,ratio::Float64,u::Float64,g::Float64)
        return new(ind,ratio,u,g)
    end
end

mutable struct Parallel_row_group
    ls::Array{Parallel_row,1}
    first::Int64
    u::Float64 # combined u
end

# find and detect all parrell constraints
# e.g., a' * x and ? * a' * x
# not only eliminate duplicate bound constraints
# algorithm: scale all vectors so first element is one and sort list
mutable struct Clever_Symmetric_KKT_solver <: abstract_KKT_system_solver
    # abstract_KKT_system_solver
    ls_solver::abstract_linear_system_solver
    factor_it::Class_iterate
    delta_x_vec::Array{Float64,1}
    delta_s_vec::Array{Float64,1}
    rhs::System_rhs
    dir::Class_point
    kkt_err_norm::Class_kkt_error
    rhs_norm::Float64
    pars::Class_parameters
    schur_diag::Array{Float64,1}
    true_x_diag::Array{Float64,1}
    diag_rescale::Array{Float64,1}

    ready::Symbol

    # Symmetric_KKT_solver only
    Q::SparseMatrixCSC{Float64,Int64}

    # these are used to eliminate parrellel constraints from the linear system solve
    para_row_info::Array{Parallel_row_group,1}
    first_para_indicies::Array{Int64,1}

    function Clever_Symmetric_KKT_solver()
      return new()
    end
end

function initialize!(kkt_solver::Clever_Symmetric_KKT_solver, intial_it::Class_iterate)
    # call this before running using the kkt_solver
    initialize!(kkt_solver.ls_solver)
    kkt_solver.dir = zero_point(dim(intial_it),ncon(intial_it))

    #start_advanced_timer(timer, "symmetric/form_system/compute_indicies");
    kkt_solver.first_para_indicies, kkt_solver.para_row_info = compute_indicies(get_jac(intial_it))
    #pause_advanced_timer(timer, "symmetric/form_system/compute_indicies");
end

function columns_are_same(A::SparseMatrixCSC{Float64,Int64},i::Int64,j::Int64)
    a_i = A[:,i]
    a_j = A[:,j]
    if a_i.nzind == a_j.nzind
        if length(a_i.nzval) == length(a_j.nzval)
            if length(a_i.nzval) > 0
                ratio = a_i.nzval[1] / a_j.nzval[1]
                if LinearAlgebra.norm(a_i - a_j * ratio,2) < 1e-16 #&& LinearAlgebra.norm(a_i,2) > 1e-16 && LinearAlgebra.norm(a_j,2) > 1e-16
                    return true
                else
                    return false
                end
            else
                println("zero elements in column")
                warn("zero elements in column")
                return true
            end
        else
            return false
        end
    else
        return false
    end
end

function rescale_cols(A::SparseMatrixCSC{Float64,Int64})
    A_rescaled = deepcopy(A)
    new_vals = A_rescaled.nzval
    i = 1
    for col = 1:size(A,2)
        a_col = A[:,col]
        if length(a_col.nzind) > 0
            tmp = a_col.nzval / a_col.nzval[1]
            l = length(tmp)
            new_vals[i:(i+l-1)] = tmp
            i += l
        end
    end

    return A_rescaled
end

function compare_columns(A::SparseMatrixCSC{Float64,Int64},i::Int64,j::Int64)
    a_i_nzind = A[:,i].nzind
    a_j_nzind = A[:,j].nzind
    if length(a_j_nzind) == 0
        return false
    elseif length(a_i_nzind) == 0
        return true
    end

    a_i_first_ind = a_i_nzind[1]
    a_j_first_ind = a_j_nzind[1]

    if a_i_first_ind != a_j_first_ind
        if a_j_first_ind < a_i_first_ind
            return true
        else
            return false
        end
    else
        if length(a_i_nzind) < length(a_j_nzind)
            return true
        elseif length(a_i_nzind) > length(a_j_nzind)
            return false
        end

        for w = 1:length(a_i_nzind)
            if a_i_nzind[w] < a_j_nzind[w]
                return true
            elseif a_i_nzind[w] > a_j_nzind[w]
                return false
            end
        end

        a_i = A[:,i]
        a_j = A[:,j]
        for k in a_i.nzind
            if a_i[k] < a_j[k]
                return true
            elseif a_i[k] > a_j[k]
                return false
            end
        end
        if i < j
            return true
        else
            return false
        end
    end
end

function sorted_col_list(A::SparseMatrixCSC{Float64,Int64})
    m = size(A,2)
    rescaled_A = rescale_cols(A)
    function comp(i::Int64,j::Int64)
        return compare_columns(rescaled_A,i,j)
    end
    sorted_cols = collect(1:m)
    sort!(sorted_cols,lt=comp)
    return sorted_cols
end

function compute_breakpoints(A::Adjoint{Float64,SparseMatrixCSC{Float64,Int64}},sorted_cols::Array{Int64,1})
    break_points = Array{Int64,1}()
    for bp = 1:length(sorted_cols)
        if bp == 1
            push!(break_points,bp)
        else
            prev = sorted_cols[bp-1]
            cur = sorted_cols[bp]
            if !columns_are_same(A,prev,cur)
                push!(break_points,bp)
            end
        end
    end
    return break_points
end

function compute_breakpoints(A::SparseMatrixCSC{Float64,Int64},sorted_cols::Array{Int64,1})
    break_points = Array{Int64,1}()
    for bp = 1:length(sorted_cols)
        if bp == 1
            push!(break_points,bp)
        else
            prev = sorted_cols[bp-1]
            cur = sorted_cols[bp]
            if !columns_are_same(A,prev,cur)
                push!(break_points,bp)
            end
        end
    end
    return break_points
end


function compute_indicies(J::SparseMatrixCSC{Float64,Int64})
    dic = Dict{Any,Array{Int64,1}}()
    J_T = SparseArrays.sparse(J')
    # REWRITE FASTER BY EMPLOYING SORTING FUNCTION.
    # a < b if
    # ahat = a / a[a_first_index], bhat = b / b[b_first_index]
    # find first i with ahat[i] < bhat[i]
    sorted_cols = sorted_col_list(J_T)

    break_points = compute_breakpoints(J_T,sorted_cols)

    no_para_indicies = sort(sorted_cols[break_points])
    para_row_info = Array{Parallel_row_group,1}()
    for k = 1:length(break_points)
        bp = break_points[k]
        if k < length(break_points)
            end_bp = break_points[k+1]
        else
            end_bp = length(sorted_cols)+1
        end
        first = sorted_cols[bp]
        ind_ls = sorted_cols[bp:(end_bp-1)]
        ls = Array{Parallel_row,1}()
        for i in ind_ls
            if i == ind_ls[1]
                ratio = 1.0
            else
                if length(J_T[:,i].nzval) > 0
                    ratio = J_T[:,i].nzval[1] / J_T[:,ind_ls[1]].nzval[1]
                else
                    ratio = 1.0
                end
                if ratio == 0.0 || isnan(ratio) || isinf(ratio)
                    println("ERROR in clever_symmetric.jl: ratio = $ratio")
                    @show ind_ls[1], i
                    @show J_T[:,i].nzval[1]
                    @show J_T[:,ind_ls[1]].nzval[1]
                    error("clever_symmetric.jl: ratio = $ratio")
                end
            end
            push!(ls, Parallel_row(i,ratio,NaN,NaN))
        end

        prg = Parallel_row_group(ls,first,NaN)
        push!(para_row_info, prg)
    end

    sort!(para_row_info,by=x->x.first)
    return no_para_indicies, para_row_info
end

function check_para_rows()
    # TODO
end

function compute_indicies(J::SparseMatrixCSC{Float64,Int64}, diag_vals::Array{Float64,1})
    ##################################################################################################################
    # INPUT
    # J::SparseMatrixCSC{Float64,Int64}
    # diag_vals::Array{Float64,1}
    #
    # OUTPUT
    #  no_para_indicies::Array{Int64,1} is a sorted list of rows that are not parallel to any other rows
    #  para_row_info::Array{} is a list containing each group of rows that are parallel.
    ##################################################################################################################

    no_para_indicies, para_row_info = compute_indicies(J)
    update_indicies!(para_row_info,diag_vals)
    return no_para_indicies, para_row_info
end

function update_indicies!(para_row_info::Array{Parallel_row_group,1}, diag_vals::Array{Float64,1})
    for row_group in para_row_info
        ls = row_group.ls
        u_inv = 0.0
        for row in ls
            row.u = diag_vals[row.ind]
            u_inv += (row.ratio)^2 * (row.u)^(-1.0)
        end
        if !(u_inv > 0.0)
            error("clever_symmetric.jl: u_inf = $u_inv !> 0.0")
        end

        group_u = 1.0/u_inv
        if !(group_u > 0.0) && !isinf(group_u)
            error("clever_symmetric.jl: $group_u = group_u !> 0.0")
        end
        row_group.u = group_u

        for row in ls
            row.g = group_u * row.ratio * (row.u)^(-1.0)
        end
    end
end

#=
A = J[lb,:]
B = J[ineq,:]
D_l, D_u, D_eq

D_bd = (D_1.^(-1) + D_2.^(-1))^(-1)

H A' B'
A D_bd 0
B 0 D_eq
=#

###########################
# code to rescale system
###########################
function create_diag_rescale_identity(factor_it::Class_iterate,u_new::Vector)
    return ones(length(factor_it.point.x) + length(u_new))
end

function create_diag_rescale_u_new(factor_it::Class_iterate,u_new::Vector)
    return [ones(length(factor_it.point.x)); factor_it.point.mu ./ sqrt.(u_new)];
end

function create_diag_rescale_u_new_dir(factor_it::Class_iterate,u_new::Vector,dir::Class_point)
    return [(LinearAlgebra.norm(dir.x,Inf) + 1e-8) * ones(length(factor_it.point.x)); factor_it.point.mu ./ sqrt.(u_new)];
end

function create_diag_rescale_u_new_and_x(factor_it::Class_iterate,u_new::Vector)
    x = factor_it.point.x
    return [ones(length(x)) / sqrt(1.0 + LinearAlgebra.norm(x,Inf)); factor_it.point.mu ./ sqrt.(u_new)];
end

function apply_rescale_to_matrix(diag_rescale,Q)
    #return spdiagm(diag_rescale) * Q * spdiagm(diag_rescale)
    return sparse(Diagonal(diag_rescale)) * Q * sparse(Diagonal(diag_rescale))
end

function apply_rescale_to_rhs(diag_rescale::Vector,clever_sym_rhs::Vector)
    return clever_sym_rhs .* diag_rescale
end

function unscale_directions(diag_rescale::Vector,dir_x_and_y::Vector)
    return dir_x_and_y .* diag_rescale
end


function form_system!(kkt_solver::Clever_Symmetric_KKT_solver, iter::Class_iterate, timer::class_advanced_timer)
    #iter.point.x = [0.7172101662159924]
    #iter.point.x = [0.2780016069742466]
    #iter.point.x = [0.0]
    #iter.point.y = [1.4998500149985006]
    #iter.point.s = [11.489889426994338]

    k = kkt_solver
    start_advanced_timer(timer, "symmetric/form_system");
    J = get_jac(iter) #.cache.J
    u = iter.point.s ./ iter.point.y;
    start_advanced_timer(timer, "symmetric/form_system/update_indicies");
    #first_para_indicies, para_row_info = compute_indicies(J,u) # ONLY DO THIS ONCE???
    #kkt_solver.first_para_indicies = first_para_indicies
    first_para_indicies = kkt_solver.first_para_indicies
    para_row_info = kkt_solver.para_row_info
    @assert( all(u.> 0.0) )
    update_indicies!(para_row_info, u)

    pause_advanced_timer(timer, "symmetric/form_system/update_indicies");


    J_new = J[first_para_indicies,:]
    #J_T_new = J_new'
    m_new = length(first_para_indicies)
    u_new = zeros(m_new) # the new diagonal elements
    for i = 1:m_new
        u_new[i] = para_row_info[i].u
    end

    if !all(u_new .>= 0)
        error("clever_symmetric.jl: minimum(u_new) = $(minimum(u_new)) < 0.0")
    end

    start_advanced_timer(timer, "symmetric/form_system/M");
    #U_new = spdiagm(u_new)
    U_new = sparse(Diagonal(u_new))
    n_var = size(J_new,2)
    M = [[get_lag_hess(iter) spzeros(n_var,m_new)]
        [J_new -U_new]];
    pause_advanced_timer(timer, "symmetric/form_system/M");

    pause_advanced_timer(timer, "symmetric/form_system");

    start_advanced_timer(timer, "symmetric/allocate");
    kkt_solver.Q = M;
    kkt_solver.factor_it = iter;

    # rescaling of matrix
    kkt_system_rescale = kkt_solver.pars.kkt.kkt_system_rescale
    if kkt_system_rescale == :none
        kkt_solver.diag_rescale = create_diag_rescale_identity(kkt_solver.factor_it,u_new)
    elseif kkt_system_rescale == :u_only
        kkt_solver.diag_rescale = create_diag_rescale_u_new(kkt_solver.factor_it,u_new)
        #kkt_solver.diag_rescale = create_diag_rescale_u_new_dir(kkt_solver.factor_it,u_new,kkt_solver.dir)
    elseif kkt_system_rescale == :u_and_x
        kkt_solver.diag_rescale = create_diag_rescale_u_new_and_x(kkt_solver.factor_it,u_new)
    end
    kkt_solver.Q = apply_rescale_to_matrix(kkt_solver.diag_rescale,kkt_solver.Q)

    kkt_solver.schur_diag = compute_schur_diag(iter)
    kkt_solver.true_x_diag = diag(M)[1:dim(iter)]
    kkt_solver.ready = :system_formed
    #println("#############################################kkt_solver.Q: ", Matrix(kkt_solver.Q))
    #println("#####################################kkt_solver.factor_it: ", kkt_solver.factor_it)
    #println("####################################kkt_solver.schur_diag: ", kkt_solver.schur_diag)

    #println("#############################################iter.point.x: ", iter.point.x)
    #println("#####################################iter.point.y: ", iter.point.y)
    #println("####################################iter.point.s: ", iter.point.s)

    pause_advanced_timer(timer, "symmetric/allocate");

end

function factor_implementation!(kkt_solver::Clever_Symmetric_KKT_solver, timer::class_advanced_timer)
    num_var = dim(kkt_solver.factor_it)
    num_con_reduced = length(kkt_solver.para_row_info)
    @assert(size(kkt_solver.Q,2) == num_var + num_con_reduced)
    return ls_factor!(kkt_solver.ls_solver, kkt_solver.Q, num_var, num_con_reduced, timer)
end

function ls_solve(mat, solver, my_rhs::Array{Float64,1}, timer::class_advanced_timer, num_ref::Int64)
	# do iterative_refinement for num_ref iterations
	sol = zeros(length(my_rhs))
	for i = 1:num_ref
		if i == 1
			err = deepcopy(my_rhs)
		else
			err = my_rhs - vector_product(mat,sol)
		end
		sol .+= ls_solve(solver, err, timer)
		#@show i, LinearAlgebra.norm(my_rhs - mat * sol)
	end
	return sol
end

function compute_direction_implementation!(kkt_solver::Clever_Symmetric_KKT_solver, timer::class_advanced_timer)
    start_advanced_timer(timer, "symmetric")
    y_org = get_y(kkt_solver.factor_it)
    s_org = get_s(kkt_solver.factor_it)
    rhs = kkt_solver.rhs;
    #println("+++++++++++++++++++++++++++++++++++y_org: ", y_org)
    #println("+++++++++++++++++++++++++++++++++++s_org: ", s_org)

    start_advanced_timer(timer, "symmetric/rhs_clever")
    symmetric_primal_rhs = rhs.primal_r + rhs.comp_r ./ y_org
    #println("+++++++++++++++++++++++++++++++++++rhs.primal_r: ", rhs.primal_r)
    #println("+++++++++++++++++++++++++++++++++++rhs.dual_r: ", rhs.dual_r)
    #println("+++++++++++++++++++++++++++++++++++rhs.comp_r: ", rhs.comp_r)
    #println("+++++++++++++++++++++++++symmetric_primal_rhs: ", symmetric_primal_rhs)

    m = length(kkt_solver.para_row_info)
    clever_symmetric_primal_rhs = zeros(m)
    for i = 1:m
        for row in kkt_solver.para_row_info[i].ls
            j = row.ind
            clever_symmetric_primal_rhs[i] += row.g * symmetric_primal_rhs[j]
        end
    end

    clever_symmetric_rhs = [rhs.dual_r; clever_symmetric_primal_rhs];
    rescaled_clever_symmetric_rhs = apply_rescale_to_rhs(kkt_solver.diag_rescale,clever_symmetric_rhs);
    pause_advanced_timer(timer, "symmetric/rhs_clever")

    # DO ITERATIVE REFINEMENT!
    #@show kkt_solver.ready

    if false
        scaled_dir_x_and_y = ls_solve(kkt_solver.ls_solver, rescaled_clever_symmetric_rhs, timer)
        #display(full(kkt_solver.Q))
        #@show LinearAlgebra.norm(clever_symmetric_rhs - vector_product(kkt_solver.Q,dir_x_and_y))
        #@show LinearAlgebra.norm(clever_symmetric_rhs - kkt_solver.Q * dir_x_and_y)
    else
        num_ref = kkt_solver.pars.kkt.ItRefine_Num # number of iterative refinement iterations
        scaled_dir_x_and_y = ls_solve(kkt_solver.Q, kkt_solver.ls_solver, rescaled_clever_symmetric_rhs, timer, num_ref)
    end
    dir_x_and_y = unscale_directions(kkt_solver.diag_rescale,scaled_dir_x_and_y);

    dir = kkt_solver.dir;
    dir.x = dir_x_and_y[1:length(rhs.dual_r)];
    v = dir_x_and_y[(1+length(rhs.dual_r)):end];

    start_advanced_timer(timer, "symmetric/dir_clever")
    y = zeros(length(dir.y))
    y_diag_vals = s_org ./ y_org;
    for j = 1:length(y)
        y[j] = y_diag_vals[j]^(-1.0) * symmetric_primal_rhs[j]
    end

    for i = 1:m
        row_group = kkt_solver.para_row_info[i]
        tmp = -(clever_symmetric_primal_rhs[i] + row_group.u * v[i])

        for row in row_group.ls
            j = row.ind
            y[j] += (row.u)^(-1.0) * row.ratio * tmp
        end
    end
    dir.y = y
    # dir.y = ...
    # dir.s = ( rhs.comp_r - dir.y .* s_org ) ./ y_org
    dir.s = eval_jac_prod(kkt_solver.factor_it, dir.x) - rhs.primal_r

    check_for_nan(dir)
    pause_advanced_timer(timer, "symmetric/dir_clever")

    start_advanced_timer(timer, "symmetric/kkt_err");
    update_kkt_error!(kkt_solver, Inf, timer)
    pause_advanced_timer(timer, "symmetric/kkt_err");
    pause_advanced_timer(timer, "symmetric")
end

function update_delta_vecs!(kkt_solver::Clever_Symmetric_KKT_solver, delta_x_vec::Array{Float64,1}, delta_s_vec::Array{Float64,1}, timer::class_advanced_timer)
    start_advanced_timer(timer, "symmetric")
    start_advanced_timer(timer, "symmetric/delta_vecs")
    kkt_solver.delta_x_vec = delta_x_vec
    kkt_solver.delta_s_vec = delta_s_vec

    if sum(abs.(delta_s_vec)) > 0.0
        error("not implemented")
    else
        if sum(abs.(delta_x_vec)) > 0.0
            n = length(delta_x_vec)
            m = size(kkt_solver.Q,2) - n
            #kkt_solver.Q = kkt_solver.Q + spdiagm([rand(length(delta_x_vec)); zeros(m)])
            kkt_solver.Q = kkt_solver.Q + sparse(Diagonal([rand(length(delta_x_vec)); zeros(m)]))
            for i = 1:length(delta_x_vec)
                kkt_solver.Q[i,i] = kkt_solver.true_x_diag[i] + delta_x_vec[i]
            end
        end

    end

    kkt_solver.ready = :delta_updated
    pause_advanced_timer(timer, "symmetric/delta_vecs")
    pause_advanced_timer(timer, "symmetric")
end
