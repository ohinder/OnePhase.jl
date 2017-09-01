#function one_phase(intial_point::Class_point)
#    return K_matrix(intial_point)
#end
#cholfact(spzeros(10,10) - speye(10))

function toy_LP()
      # min 10 x, x >= 0
    intial_point = Class_point([1.0], [1.0], [1.0], 1.0)

    A = ones(1,1)
    b = [0.0]
    c = [10.0]

    LP = Construct_LP_problem(A, b, c);

    intial_it = Class_iterate(intial_point, LP, Class_local_info(0.0));

    return intial_it
end

function toy_LP2()
    # min 4 x + y,
    # x + y >= 0
    # 5 * x + y >= 0
    # x + 5 * y >= 0

    x = [1.0; 1.0]

    A = [[1.0 1.0]; [5.0 1.0]; [1.0 5.0]]
    b = [0.0; 0.0; 0.0];
    c = [4.0; 1.0];
    LP = Construct_LP_problem(A, b, c);

    s = A * x - b + 2.0
    mu = 1.0
    y = mu ./ s

    intial_point = Class_point(x, y, s, mu)

    intial_it = Class_iterate(intial_point, LP, Class_local_info(0.0));

    return intial_it
end
