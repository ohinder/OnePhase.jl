type interval
  
end

function bisection(neg_f_position::Float64, pos_f_position::Float64, func::Function, eps1::Float64)
    return bisection(neg_f_position, pos_f_position, func, eps1, 0.0)
end

function bisection(neg_f_position::Float64, pos_f_position::Float64, func::Function, eps1::Float64, eps2::Float64)
    func_neg_f_position = func(neg_f_position)
    if (func_neg_f_position >= 0.0)
        @show func_neg_f_position
        error("func(neg_f_position) >= 0.0")
    end

    func_pos_f_position = func(pos_f_position)
    if (func_pos_f_position <= 0.0)
        @show func_pos_f_position
        @show func(pos_f_position)
        error("func(pos_f_position) <= 0.0")
    end

    return recursive_bisection(neg_f_position, pos_f_position, func, eps1, eps2)
end


function recursive_bisection(neg_f_position::Float64, pos_f_position::Float64, func::Function, eps1::Float64, eps2::Float64)
    # find zero of increasing func
    # do with for loop?
    centre_point = (pos_f_position + neg_f_position)/2.0
    val = func(centre_point)

    terminate = abs(val) > eps1 && abs(neg_f_position - pos_f_position) > eps2

    if val > 0.0
        pos_f_position = centre_point
    elseif val < 0.0
        neg_f_position = centre_point
    elseif val == 0.0
        return centre_point, centre_point
    end


    if abs(val) < eps1
        return centre_point, centre_point
    elseif abs(neg_f_position - pos_f_position) < eps2
        return neg_f_position, pos_f_position
    else
        return recursive_bisection(neg_f_position, pos_f_position, func, eps1, eps2)
    end
end

if false
    function my_test_func(x::Float64)
        return x
    end

    bisection(-1.0,2.0,my_test_func,0.01)
end
