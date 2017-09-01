#using MAT

function if_mkdir(dir::String)
  if !isdir(dir)
     mkdir(dir)
  end
end

function rd(num::Float64, digits::Int64=2) # round
    if digits == 2
      return pd(@sprintf("%.2e", num))
    elseif digits == 3
      return pd(@sprintf("%.3e", num))
    elseif digits == 4
      return pd(@sprintf("%.4e", num))
    elseif digits == 5
      return pd(@sprintf("%.5e", num))
    else
      error("digits $digits not supported for rd")
    end
end

function pd(input)
    return pd(string(input))
end

function pd(str::String)
    return rpad(str, 11 ," ")
end

function pd(input, pad_len::Int64)
    return pd(string(input), pad_len)
end

function pd(str::String, pad_len::Int64)
    return rpad(str, pad_len ," ")
end

function solve_quadratic(a::Float64, b::Float64, c::Float64)
    val = -b/(2 * a)
    pm = sqrt(b^2 - 4 * a * c) / (2 * a)

    return (val + pm, val - pm)
end

function my_warn(str::String)
    println("WARNING: ", str)
end

function densest_row(mat) #::SparseMatrixCSC{Float64,Int32})
    density = 0
    for i = 1:size(mat,1)
      density = max(density, nnz(mat[i,:]))
    end

    return density
end

function densest_col(mat) #::SparseMatrixCSC{Float64,Int32})
    density = 0
    for i = 1:size(mat,2)
      density = max(density, nnz(mat[:,i]))
    end

    return density
end
