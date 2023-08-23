@testset "apply!" begin
    @testset "FermiCreationOperator" begin
        state = falses(5)
        op = FermiCreationOperator(3)
        apply!(state, op) == 1.0
        @test state == [false, false, true, false, false]
    end
    
    @testset "FermiAnnihilationOperator" begin
        state = [false, false, true, false, false]
        op = FermiAnnihilationOperator(3)
        apply!(state, op) == 1.0
        @test state == falses(5)
    end
    
    @testset "FermiNumberOperator" begin
        state = [false, false, true, false, false]
        op = FermiNumberOperator(3)
        @test apply!(state, op) == 1.0
        @test state == [false, false, true, false, false]
    end
    
    @testset "FermiHoppingOperator" begin
        state = [false, false, false, false, true]
        op = FermiHoppingOperator(3, 5)
        apply!(state, op) == 1.0
        @test state == [false, false, true, false, false]
    end

end


@testset "cdcd == 0" begin
    @testset "AscendingTwoLevelBasis" begin
        basis = AscendingTwoLevelBasis(1)
        cd = spmatrix(FermiCreationOperator(1), basis)
        @test iszero(cd*cd)

        for L in 2:5
            basis = AscendingTwoLevelBasis(L)
            cd = spmatrix(FermiCreationOperator(2), basis)
            @test iszero(cd*cd)
        end

        basis = AscendingTwoLevelBasis(6)
        for i in 1:6
            c = spmatrix(FermiCreationOperator(i), basis)
            @test iszero(c*c)
        end
    end

    @testset "NumberConservedBasis" begin
        for L in 2:5
            basis = NumberConservedBasis(L, ceil(Int,L/2))
            cd = spmatrix(FermiCreationOperator(2), basis)
            @test iszero(cd*cd)
        end
    end
end

@testset "cc == 0" begin
    @testset "AscendingTwoLevelBasis" begin
        basis = AscendingTwoLevelBasis(1)
        c = spmatrix(FermiAnnihilationOperator(1), basis)
        @test iszero(c*c)

        for L in 2:5
            basis = AscendingTwoLevelBasis(L)
            c = spmatrix(FermiAnnihilationOperator(2), basis)
            @test iszero(c*c)
        end

        basis = AscendingTwoLevelBasis(6)
        for i in 1:6
            c = spmatrix(FermiAnnihilationOperator(i), basis)
            @test iszero(c*c)
        end
    end

    @testset "NumberConservedBasis" begin
        for L in 2:5
            basis = NumberConservedBasis(L, ceil(Int,L/2))
            c = spmatrix(FermiAnnihilationOperator(2), basis)
            @test iszero(c*c)
        end
    end
end

@testset "c == 0 for NumberConservedBasis" begin 
    for L in 1:5
        basis = NumberConservedBasis(L, ceil(Int,L/2))
        for i in 1:L
            c = spmatrix(FermiAnnihilationOperator(i), basis)
            @test iszero(c)
        end
    end
end

@testset "cd == 0 for NumberConservedBasis" begin 
    for L in 1:5
        basis = NumberConservedBasis(L, ceil(Int,L/2))
        for i in 1:L
            cd = spmatrix(FermiCreationOperator(i), basis)
            @test iszero(cd)
        end
    end
end

@testset "cdc == n" begin
    @testset "TensoAscendingTwoLevelBasisrBasis" begin
        basis = AscendingTwoLevelBasis(1)
        cd = spmatrix(FermiCreationOperator(1), basis)
        c = spmatrix(FermiAnnihilationOperator(1), basis)
        n = spmatrix(FermiNumberOperator(1), basis)
        @test cd*c == n

        for L in 2:5
            basis = AscendingTwoLevelBasis(L)
            cd = spmatrix(FermiCreationOperator(2), basis)
            c = spmatrix(FermiAnnihilationOperator(2), basis)
            n = spmatrix(FermiNumberOperator(2), basis)
            @test cd*c == n
        end
        
    end
end

@testset "cd_i c_j == FermiHopping(i,j)" begin 
    @testset "AscendingTwoLevelBasis" begin
        for L in 2:5
            for i in 1:L
                for j in 1:L
                    i == j && continue 

                    basis = AscendingTwoLevelBasis(L)
                    cd = spmatrix(FermiCreationOperator(i), basis)
                    c = spmatrix(FermiAnnihilationOperator(j), basis)
                    h = spmatrix(FermiHoppingOperator(i,j), basis)
                    @test cd*c == h
                end
            end
        end
    end
end