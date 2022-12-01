module QuantumOperators

using Reexport
using LinearAlgebra
using SparseArrays

@reexport using QuantumBases

export spmatrix

# Spins
export SigmaX, SigmaY, SigmaZ, SigmaPlus, SigmaMinus
export SigmaXX, SigmaYY, SigmaZZ
export SigmaPlusMinus, SigmaMinusPlus

# Spinless fermions
export FermiNumberOperator, FermiHoppingOperator

include("abstract.jl")
include("spin.jl")
include("fermion.jl")

end
