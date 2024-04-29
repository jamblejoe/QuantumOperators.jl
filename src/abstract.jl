abstract type AbstractOperator{vT} end

Base.eltype(::AbstractOperator{vT}) where {vT} = vT


function spmatrix(op::AbstractOperator, basis::AbstractBasis)
    spmatrix(op, basis, basis)
end

function spmatrix(op::AbstractOperator{vT}, 
    basis1::AbstractBasis, basis2::AbstractBasis) where {vT}

    basis1.L == basis2.L || error("basis1 and basis2 must have same number of sites. Got $(basis1.L) and $(basis2.L)")
    L = basis1.L

    eltype(basis1) == eltype(basis2) || error("basis1 and basis2 must to have the same eltype!")

    #basis_element = BitVector(undef, L)
    basis_element = eltype(basis1)(undef, L)

    rows = Int[]
    cols = Int[]
    values = vT[]

    for i in eachindex(basis1)
        getstate!(basis_element, basis1, i)

        #basis_element[site_index] = !basis_element[site_index]
        val = apply!(basis_element, op)
        iszero(val) && continue

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