NORM_TYPE = Inf

function lb_s_thres(iter::Class_iterate,dir::Class_point,pars::Class_parameters)
  ex = pars.ls.fraction_to_boundary_predict_exp
  #alpha_P_tilde = simple_max_step(iter.point.s,dir.s,zeros(length(dir.s)))
  #x_thres = max.(2.0 * norm(dir.x,NORM_TYPE)^(2 + ex), norm(dir.x,NORM_TYPE)^(1.0 + ex) + norm(dir.y,NORM_TYPE)^0.5 * norm(dir.x,NORM_TYPE))
  x_thres = norm(dir.x,NORM_TYPE)^ex #+ norm(dir.y,NORM_TYPE) + get_delta(iter)

  return min.(iter.point.s, norm(dir.x,NORM_TYPE) * x_thres)
end

function lb_s_predict(iter::Class_iterate,dir::Class_point,pars::Class_parameters)
  #@show iter.frac_bd_predict
  return iter.frac_bd_predict .* lb_s_thres(iter,dir,pars)
end

function lb_y(iter::Class_iterate,dir::Class_point,pars::Class_parameters)
  #alpha_P_tilde = simple_max_step(iter.point.s,dir.s,zeros(length(dir.s)))
  return iter.frac_bd .*  iter.point.y * min(1.0, norm(dir.x,NORM_TYPE) )
end

function lb_s(iter::Class_iterate,dir::Class_point,pars::Class_parameters)
    #alpha_P_tilde = simple_max_step(iter.point.s,dir.s,zeros(length(dir.s)))
    lb_s = iter.frac_bd .* lb_s_thres(iter,dir,pars)
    #min.(iter.point.s, norm(dir.x,NORM_TYPE)^2)
    #@show lb_s
    return lb_s
end

function simple_max_step(val::Array{Float64,1}, dir::Array{Float64,1}, lb::Array{Float64,1})
    ratio = maximum( [1.0; -dir ./ (val - lb) ] ) #; -q * dir.y ./ iter.point.y] );
    max_step = 1.0 / ratio
    #@show max_step
    return max_step
end
