"""
Pauli matrices
σ_x = [0 1; 1 0]
σ_y = [0 -im; im 0]
σ_z = [1 0; 0 -1]

σ_+ = [0 1; 0 0] = 1/2 (σ_x + iσ_y)
σ_- = [0 0; 1 0] = 1/2 (σ_x - iσ_y)

|0> = [0, 1]' == spin down
|1> = [1, 0]' == spin up
"""


struct SingleBodyOperator{vT, sT} <: AbstractOperator{vT}
    site::Int
end

#=
abstract type SigmaX <: AbstractOperator end
abstract type SigmaY <: AbstractOperator end
abstract type SigmaZ <: AbstractOperator end
abstract type SigmaPlus <: AbstractOperator end
abstract type SigmaMinus <: AbstractOperator end
=#

#struct Operator{T} end
#const SigmaX = Operator{:X}()
#const SigmaY = Operator{:Y}()
#const SigmaZ = Operator{:Z}()
#const SigmaPlus = Operator{:+}()
#const SigmaMinus = Operator{:-}()
SigmaX(i::Integer) = SingleBodyOperator{Int, :X}(i)
SigmaY(i::Integer) = SingleBodyOperator{complex(Int), :Y}(i)
SigmaZ(i::Integer) = SingleBodyOperator{Int, :Z}(i)
SigmaPlus(i::Integer) = SingleBodyOperator{Int, :+}(i)
SigmaMinus(i::Integer) = SingleBodyOperator{Int, :-}(i)



function apply!(state::BitVector, op::SingleBodyOperator{vT,:X}) where {vT}
    site = op.site
    state[site] = !state[site]
    return one(vT)
end
function apply!(state::BitVector, op::SingleBodyOperator{vT,:Y}) where {vT}
    site = op.site
    s = state[site]
    state[site] = !s
    return s ? -one(vT)*im : one(vT)*im
end
function apply!(state::BitVector, op::SingleBodyOperator{vT,:Z}) where {vT}
    site = op.site
    s = state[site]
    return s ? -one(vT) : one(vT)
end

function apply!(state::BitVector, op::SingleBodyOperator{vT,:+}) where {vT}
    site = op.site
    s = state[site]
    state[site] = !state[site]
    return s ? zero(vT) : one(vT)
end
function apply!(state::BitVector, op::SingleBodyOperator{vT,:-}) where {vT}
    site = op.site
    s = state[site]
    state[site] = !state[site]
    return s ? one(vT) : zero(vT)
end

spmatrix(::SingleBodyOperator{vT,:X}) where {vT} = sparse(vT[0 1; 1 0])
spmatrix(::SingleBodyOperator{vT,:Y}) where {vT} = sparse(vT[0 -im; im 0])
spmatrix(::SingleBodyOperator{vT,:Z}) where {vT} = sparse(vT[1 0; 0 -1])
spmatrix(::SingleBodyOperator{vT,:+}) where {vT} = sparse(vT[0 1; 0 0])
spmatrix(::SingleBodyOperator{vT,:-}) where {vT} = sparse(vT[0 0; 1 0])

struct TwoBodyOperator{vT,T1,T2} <: AbstractOperator{vT}
    op1::SingleBodyOperator{vT,T1}
    op2::SingleBodyOperator{vT,T2}
end

# the order is important when i=j
# σ_i^+ σ_j^-
#SigmaPlusMinus(i::Integer, j::Integer) = TwoBodyOperator{:±}(i,j)
SigmaPlusMinus(i::Integer, j::Integer) = TwoBodyOperator(SigmaPlus(i),SigmaMinus(j))
# σ_i^- σ_j^+
#SigmaMinusPlus(i::Integer, j::Integer) = TwoBodyOperator{:∓}(i,j)
SigmaMinusPlus(i::Integer, j::Integer) = TwoBodyOperator(SigmaMinus(i),SigmaPlus(j))

SigmaXX(i::Integer, j::Integer) = TwoBodyOperator(SigmaX(i), SigmaX(j))
SigmaYY(i::Integer, j::Integer) = TwoBodyOperator(SigmaY(i), SigmaY(j))
SigmaZZ(i::Integer, j::Integer) = TwoBodyOperator(SigmaZ(i), SigmaZ(j))

# create a spin at site 1 and annihilate a spin at site 2
function apply!(state::BitVector, op::TwoBodyOperator{vT, :+, :-}) where {vT}
    site1 = op.op1.site
    site2 = op.op2.site

    if site1 == site2
        return state[site1] ? one(vT) : zero(vT)
    else
        # returns true if and only if site 1 is empty and site 2 is occupied
        s = !state[site1] && state[site2]
        state[site1] = one(vT)
        state[site2] = zero(vT)
        return s ? one(vT) : zero(vT)
    end
end

# annihilate a spin at site 1 and create a spin at site 2
function apply!(state::BitVector, op::TwoBodyOperator{vT,:-,:+}) where {vT}
    site1 = op.op1.site
    site2 = op.op2.site

    if site1 == site2
        return state[site1] ? zero(vT) : one(vT)
    else
        # returns true if and only if site 1 is occupied and site 2 is empty
        s = state[site1] && !state[site2]
        state[site1] = zero(vT)
        state[site2] = one(vT)
        return s ? one(vT) : zero(vT)
    end
end

#=
function spmatrix(op::AbstractOperator, basis::AbstractBasis, T::Type=Float64)
    spmatrix(op, basis, basis, T)
end

function spmatrix(op::AbstractOperator, 
    basis1::AbstractBasis, basis2::AbstractBasis,
    T::Type=Float64)

    basis1.L == basis2.L || error("basis1 and basis2 must have same number of sites. Got $(basis1.L) and $(basis2.L)")
    L = basis1.L

    basis_element = BitVector(undef, L)

    rows = Int[]
    cols = Int[]
    values = T[]

    for i in eachindex(basis1)
        getstate!(basis_element, basis1, i)

        #basis_element[site_index] = !basis_element[site_index]
        val = apply!(basis_element, op, T)

        if basis_element in basis2

            # get the index of the new basis element
            j = getposition(basis2, basis_element)

            # store the indices
            push!(rows, j)
            push!(cols, i)

            # calculate and store the value
            push!(values, val)

        end
    end

    # create the sparse matrices from the rows/columns/values
    spm = sparse(rows, cols, values, length(basis2), length(basis1))

    return spm
end
=#

####################################################################
#
# TensorBasis specialized functions
#
####################################################################


function spmatrix(op::SingleBodyOperator, basis::TensorBasis)
    L = basis.L
    site_index = op.site
    1 <= site_index <= L || error("site must be in [1,$L], got $i")

    m = spmatrix(op)
    M = kron(I(2^(site_index-1)), m, I(2^(L-site_index)))
    M
end


function spmatrix(op::TwoBodyOperator, basis::TensorBasis)
    L = basis.L
    op1 = op.op1
    op2 = op.op2
    site1 = op1.site
    site2 = op2.site

    1 <= site1 <= L || error("site1 must be in [1,$L], got $(site1)")
    1 <= site2 <= L || error("site2 must be in [1,$L], got $(site2)")


    m1 = spmatrix(op1)
    m2 = spmatrix(op2)

    M = nothing
    if site1 < site2
        M = kron(I(2^(site1-1)), m1, I(2^(site2-site1-1)), m2, I(2^(L-site2)))
        #kron(I(2^(L-site2)), m2, I(2^(site2-site1-1)), m1, I(2^(site1-1)))
    elseif site2 < site1
        M = kron(I(2^(site2-1)), m2, I(2^(site1-site2-1)), m1, I(2^(L-site1)))
        #kron(I(2^(L-site1)), m1, I(2^(site1-site2-1)), m2, I(2^(site2-1)))
    else
        # site1==site2
        M = kron(I(2^(site1-1)), m1*m2, I(2^(L-site1)))
        #kron(I(2^(L-site1)), m2*m1, I(2^(site1-1)))
    end

    return M
end


####################################################################
#
# Algebra
#
####################################################################

#=
Pauli = SingleBodyOperator
import Base.:*
import Base.sort!

struct SpinOperator{T<:Number}
    phase::T
    ops::Vector{SingleBodyOperator}
end
function SpinOperator(ops::Vararg{Pauli}; T::Type=ComplexF64)
    SpinOperator(one(T), [ops...])
end


*(a::Pauli, b::Pauli) = SpinOperator(a, b)


function sort!(op::SpinOperator)
    sort!(op.ops, by=x->x.site)
end


function contract(op::SpinOperator)
    sort!(op)
    for i in 1:length(op.ops)
    end
end

=#

#### impl 2

#=
Pauli = SingleBodyOperator

struct SpinOperator{T<:Number}
    phase::T
    ops::Vector{SingleBodyOperator}
end
=#