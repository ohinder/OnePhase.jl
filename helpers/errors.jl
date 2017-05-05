function my_throw(e::Exception, location_msg::String)
    println("ERROR at ", location_msg)
    throw(e)
end

#=
type exception_recursion <: Exception
    prev_e::Exception
    location_msg::String
    function exception_recursion(prev_e::Exception, location_msg::String)
        return new(prev_e, location_msg)
    end
end

function print(e::exception_recursion)
    println("ERROR at ", e.location_msg)
    print(e.prev_e)
end

function print(e::Exception)
    print(e.msg)
end
=#
