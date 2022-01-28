
"""
c_i^dagger c_i
"""
struct FermiNumberOperator <: AbstractOperator
    i::Int
end

"""
c_i^dagger c_j

"""
struct FermiHoppingOperator <: AbstractOperator
    i::Int
    j::Int
end

function apply!(state::AbstractArray{<:Bool}, op::FermiNumberOperator, T::Type=Float64)
    i = op.i
    state[i] == 0 && return zero(T)
    state[i] == 1 && return one(T)
end

function apply!(state::AbstractArray{<:Bool}, op::FermiHoppingOperator, T::Type=Float64)
    i = op.i
    j = op.j

    state[i] == 1 && return zero(T)
    state[j] == 0 && return zero(T)

    state[i] = 0
    state[j] = 1

    s = 0
    for l in (min(i,j)+1):(max(i,j)-1)
        s += state[l]
    end

    return (-one(T))^s

end

