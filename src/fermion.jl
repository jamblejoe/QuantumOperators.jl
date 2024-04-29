

"""
c_i^dagger
"""
struct FermiCreationOperator{vT} <: AbstractOperator{vT}
    i::Int
    FermiCreationOperator(i::Integer) = new{Int}(i)
end

function apply!(state::AbstractArray{<:Bool}, op::FermiCreationOperator{vT}) where {vT}
    i = op.i
    isone(state[i]) && return zero(vT)
    state[i] = 1
    
    s = 0
    for l in 1:(i-1)
        s += state[l]
    end
    return iseven(s) ? one(vT) : -one(vT)
end

"""
c_i
"""
struct FermiAnnihilationOperator{vT} <: AbstractOperator{vT}
    i::Int
    FermiAnnihilationOperator(i::Integer) = new{Int}(i)
end

function apply!(state::AbstractArray{<:Bool}, op::FermiAnnihilationOperator{vT}) where {vT}
    i = op.i
    iszero(state[i]) && return zero(vT)
    state[i] = 0
    
    s = 0
    for l in 1:(i-1)
        s += state[l]
    end
    return iseven(s) ? one(vT) : -one(vT)
end


"""
c_i^dagger c_i
"""
struct FermiNumberOperator{vT} <: AbstractOperator{vT}
    i::Int
    FermiNumberOperator(i::Integer) = new{Int}(i)
end

function apply!(state::AbstractArray{<:Bool}, op::FermiNumberOperator{vT}) where {vT}
    i = op.i
    iszero(state[i]) && return zero(vT)
    isone(state[i])  && return one(vT)
end

"""
c_i^dagger c_j

"""
struct FermiHoppingOperator{vT} <: AbstractOperator{vT}
    i::Int
    j::Int
    function FermiHoppingOperator(i::Integer, j::Integer)
        i == j && throw(DomainError("Only different sites supported. Use FermiNumberOperator instead."))
        new{Int}(i,j)
    end
end

function apply!(state::AbstractArray{<:Bool}, op::FermiHoppingOperator{vT}) where {vT}
    i = op.i
    j = op.j

    isone(state[i])  && return zero(vT)
    iszero(state[j]) && return zero(vT)

    state[i] = 1
    state[j] = 0

    s = 0
    for l in (min(i,j)+1):(max(i,j)-1)
        s += state[l]
    end

    return iseven(s) ? one(vT) : -one(vT)
end


"""
    Every computational basis state is an eigenstate of the parity operator.
    Depending on the number of particles, the parity operator has eigenvalue
    +1 or -1. It is +1 if the number of particles is even and -1 if the number
    of particles is odd.

"""
struct ParityOperator{vT} <: AbstractOperator{vT} 
    ParityOperator() = new{Int}()
end

function QuantumOperators.apply!(state::AbstractArray{<:Bool}, op::ParityOperator{vT}) where {vT}
	s = sum(state)
	iseven(s) && return  one(vT)
	isodd(s)  && return -one(vT)
end









