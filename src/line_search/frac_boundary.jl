function lb_s_predict(iter,dir,pars)
  ex = pars.ls.fraction_to_boundary_predict_exp
  x_thres = max(2.0 * norm(dir.x,Inf)^2, norm(dir.x,Inf)^ex)
  #bd = ones(length(iter.point.y)) * pars.ls.fraction_to_boundary_predict
  lb_s_predict = iter.frac_bd_predict .*  min.(iter.point.s, x_thres)
  #@show lb_s_predict
  return lb_s_predict
end

function lb_y(iter,dir,pars)
  return iter.frac_bd .*  iter.point.y * min(1, norm(dir.x,Inf) )
end

function lb_s(iter,dir,pars)
    lb_s = iter.frac_bd .*  min.(iter.point.s, norm(dir.x,Inf)^2)
    #@show lb_s
    return lb_s
end

function simple_max_step(val::Array{Float64,1}, dir::Array{Float64,1}, lb::Array{Float64,1})
    ratio = maximum( [1.0; -dir ./ (val - lb) ] ) #; -q * dir.y ./ iter.point.y] );
    max_step = 1.0 / ratio
    #@show max_step
    return max_step
end
