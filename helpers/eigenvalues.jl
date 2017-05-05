function eig_min!(linear_solver::abstract_linear_system_solver, K::SparseMatrixCSC{Float64,Int64}, x_0::Array{Float64,1}, n::Int64, m::Int64, tol::Float64)
	# computes smallest magnitude eigenvalue and eigenvector viz inverse iteration
	try
		@assert(n + m == size(K,1))
		@assert(size(K,1) == size(K,2))
		@assert(n == length(x_0))

		eigen_value = 0.0;
		x = x_0;
		sol = zeros(n + m);
		MAX_IT = 30;
		err = Inf;
		for i = 1:MAX_IT

			rhs = [x; zeros(m)];

			sol = ls_solve(linear_solver, rhs);
			x_old = x;
			x = sol[1:n];
			eigen_value = 1.0/norm(x,2)
			x = x * eigen_value;

			err = norm(x - x_old, 1)

			if err < tol
				println(eigen_value)
				println("converged!")
				break
			end
		end

		y = sol[(n+1):(n+m)];

		return eigen_value, x, y, err
	catch e
		println("ERROR in eig_min")
		throw(e)
	end
end

function is_B_matrix(B::SparseMatrixCSC{Float64,Int64})
	# checks that matrix is made of 1x1 or 2x2 blocks

end

function eigsigns_of_B(B::SparseMatrixCSC{Float64,Int64}, toler::Float64)
	#

	@assert(size(B,1) == size(B,2))
	n = size(B,1);
	trans_eigenvals = zeros(n);
	i = 1;
	while i <= n
		if i == n || B[i, i+1] == 0
			trans_eigenvals[i] = B[i, i];
			i = i + 1;
		else
			trans_eigenvals[i:(i+1)] = eigvals(full(B[i:i+1,i:i+1]));
			i = i + 2;
		end
	end

	pos_eigs = sum(trans_eigenvals .> toler);
	neg_eigs = sum(trans_eigenvals .< toler);
	zero_eigs = length(trans_eigenvals) - pos_eigs - neg_eigs;

	return pos_eigs, neg_eigs, zero_eigs
end


if false
	H = spdiagm([1.1,-1.0]);
	A = sparse([[1.0, 1.0] ])
	K = [ [H A]; [A' spzeros(1,1)] ];


	linear_solver = linear_solver_JULIA();
	linear_solver.initialize(K);
	linear_solver.ls_factor();

	x_0 = randn(2);

	@show eig_min!(linear_solver, K, x_0, 2, 1, 1e-6)
end
