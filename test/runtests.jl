using Test
using QuantumOperators
using LinearAlgebra
using SparseArrays

# these functions are not exported by QuantumOperators
# therefore we have to import them explicitly
apply! = QuantumOperators.apply!

include("spintests.jl")
include("fermiontests.jl")