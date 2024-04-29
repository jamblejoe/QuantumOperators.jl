@testset "Spin operators" begin

@testset "Pauli matrices" begin
    @test spmatrix(SigmaX(1)) == [0 1; 1 0]
    @test spmatrix(SigmaY(1)) == [0 -im; im 0]
    @test spmatrix(SigmaZ(1)) == [1 0; 0 -1]
    @test spmatrix(SigmaPlus(1)) == [0 1; 0 0]
    @test spmatrix(SigmaMinus(1)) == [0 0; 1 0]

    @test spmatrix(SigmaPlus(1)) == 1/2 * (spmatrix(SigmaX(1)) + im * spmatrix(SigmaY(1)))
    @test spmatrix(SigmaMinus(1)) == 1/2 * (spmatrix(SigmaX(1)) - im * spmatrix(SigmaY(1)))

    b = TensorBasis(1)
    @test spmatrix(SigmaX(1)) == spmatrix(SigmaX(1), b)
    @test spmatrix(SigmaY(1)) == spmatrix(SigmaY(1), b)
    @test spmatrix(SigmaZ(1)) == spmatrix(SigmaZ(1), b)
    @test spmatrix(SigmaPlus(1)) == spmatrix(SigmaPlus(1), b)
    @test spmatrix(SigmaMinus(1)) == spmatrix(SigmaMinus(1), b)

    @test spmatrix(SigmaPlus(1)) == 1/2 * (spmatrix(SigmaX(1), b) + im * spmatrix(SigmaY(1), b))
    @test spmatrix(SigmaMinus(1)) == 1/2 * (spmatrix(SigmaX(1), b) - im * spmatrix(SigmaY(1), b))

end

@testset "SigmaPlusMinus == SigmaPlus * SigmaMinus" begin
    @testset "TensorBasis" begin
        
        basis =  TensorBasis(4)
        for i in 1:3
            M1 = spmatrix(SigmaPlusMinus(i,i+1), basis)
            M2 = spmatrix(SigmaPlus(i), basis) * spmatrix(SigmaMinus(i+1), basis)
            @test M1 == M2
        end
    end

    @testset "AscendingTwoLevelBasis" begin
        
        basis =  AscendingTwoLevelBasis(4)
        for i in 1:3
            M1 = spmatrix(SigmaPlusMinus(i,i+1), basis)
            M2 = spmatrix(SigmaPlus(i), basis) * spmatrix(SigmaMinus(i+1), basis)
            @test M1 == M2
        end
    end
end


@testset "SigmaPlusMinus" begin

    @testset "NumberConservedBasis" begin
        
        M = spmatrix(SigmaPlusMinus(1,2), NumberConservedBasis(4,2))
        M_corr = spzeros(Int, 6,6)
        M_corr[4,2] = 1
        M_corr[5,3] = 1
        @test M == M_corr
        
    end

    @testset "AscendingTwoLevelBasis" begin
        
        M = spmatrix(SigmaPlusMinus(1,2), AscendingTwoLevelBasis(2))
        M_corr = spzeros(Int, 4,4)
        M_corr[3,2] = 1
        @test M == M_corr

        M = spmatrix(SigmaPlusMinus(2,1), AscendingTwoLevelBasis(2))
        M_corr = spzeros(Int, 4,4)
        M_corr[2,3] = 1
        @test M == M_corr

        M = spmatrix(SigmaPlusMinus(1,2), AscendingTwoLevelBasis(3))
        M_corr = spzeros(Int, 8,8)
        M_corr[5,3] = 1
        M_corr[6,4] = 1
        @test M == M_corr

        M = spmatrix(SigmaPlusMinus(2,1), AscendingTwoLevelBasis(3))
        M_corr = spzeros(Int, 8,8)
        M_corr[3,5] = 1
        M_corr[4,6] = 1
        @test M == M_corr

        M = spmatrix(SigmaPlusMinus(1,3), AscendingTwoLevelBasis(3))
        M_corr = spzeros(Int, 8,8)
        M_corr[5,2] = 1
        M_corr[7,4] = 1
        @test M == M_corr

        M = spmatrix(SigmaPlusMinus(3,1), AscendingTwoLevelBasis(3))
        M_corr = spzeros(Int, 8,8)
        M_corr[2,5] = 1
        M_corr[4,7] = 1
        @test M == M_corr

        M = spmatrix(SigmaPlusMinus(2,3), AscendingTwoLevelBasis(3))
        M_corr = spzeros(Int, 8,8)
        M_corr[3,2] = 1
        M_corr[7,6] = 1
        @test M == M_corr

        M = spmatrix(SigmaPlusMinus(3,2), AscendingTwoLevelBasis(3))
        M_corr = spzeros(Int, 8,8)
        M_corr[2,3] = 1
        M_corr[6,7] = 1
        @test M == M_corr
        
    end

    
end


@testset "Adjoint of operators" begin

    @testset "SigmaPlusMinus" begin
        @testset "TensorBasis($L)" for L in [3,4,5]
            b = AscendingTwoLevelBasis(L)
            @testset "i=$i, j=$j" for i in 1:b.L, j in 1:b.L
                @test spmatrix(SigmaPlusMinus(i,j),b) == adjoint(spmatrix(SigmaPlusMinus(j,i),b))
            end
        end

        @testset "AscendingTwoLevelBasis($L)" for L in [3,4,5]
            b = AscendingTwoLevelBasis(L)
            @testset "i=$i, j=$j" for i in 1:b.L, j in 1:b.L
                @test spmatrix(SigmaPlusMinus(i,j),b) == adjoint(spmatrix(SigmaPlusMinus(j,i),b))
            end
        end

        @testset "NumberConservedBasis($L, $N)" for (L,N) in [(3,1), (4,2), (6,4)]
            b = NumberConservedBasis(L,N)
            @testset "i=$i, j=$j" for i in 1:b.L, j in 1:b.L
                @test spmatrix(SigmaPlusMinus(i,j),b) == adjoint(spmatrix(SigmaPlusMinus(j,i),b))
            end
        end

    end

    @testset "SigmaMinusPlus" begin

        @testset "AscendingTwoLevelBasis($L)" for L in [3,4,5]
            b = AscendingTwoLevelBasis(L)
            @testset "i=$i, j=$j" for i in 1:b.L, j in 1:b.L
                @test spmatrix(SigmaMinusPlus(i,j),b) == adjoint(spmatrix(SigmaMinusPlus(j,i),b))
            end
        end
        
        @testset "NumberConservedBasis($L, $N)" for (L,N) in [(3,1), (4,2), (6,4)]
            b = NumberConservedBasis(L,N)
            @testset "i=$i, j=$j" for i in 1:b.L, j in 1:b.L
                @test spmatrix(SigmaMinusPlus(i,j),b) == adjoint(spmatrix(SigmaMinusPlus(j,i),b))
            end
        end
    end
    
end

end