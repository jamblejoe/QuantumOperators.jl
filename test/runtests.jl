using Test
using QuantumOperators
using SparseArrays

apply! = QuantumOperators.apply!

include("spintests.jl")
include("fermiontests.jl")