@testset "SigmaPlusMinus" begin

    @testset "NumberConservedBasis" begin
        
        M = spmatrix(SigmaPlusMinus(1,2), NumberConservedBasis(4,2), Int)
        M_corr = spzeros(Int, 6,6)
        M_corr[4,2] = 1
        M_corr[5,3] = 1
        @test M == M_corr
        
    end
    
end


@testset "Adjoint of operators" begin

    @testset "SigmaPlusMinus" begin
        @testset "TensorBasis($L)" for L in [3,4,5]
            b = TensorBasis(L)
            @testset "i=$i, j=$j" for i in 1:b.L, j in 1:b.L
                @test spmatrix(SigmaPlusMinus(i,j),b, Int) == adjoint(spmatrix(SigmaPlusMinus(j,i),b, Int))
            end
        end

        @testset "NumberConservedBasis($L, $N)" for (L,N) in [(3,1), (4,2), (6,4)]
            b = NumberConservedBasis(L,N)
            @testset "i=$i, j=$j" for i in 1:b.L, j in 1:b.L
                @test spmatrix(SigmaPlusMinus(i,j),b, Int) == adjoint(spmatrix(SigmaPlusMinus(j,i),b, Int))
            end
        end

    end

    @testset "SigmaMinusPlus" begin

        @testset "TensorBasis($L)" for L in [3,4,5]
            b = TensorBasis(L)
            @testset "i=$i, j=$j" for i in 1:b.L, j in 1:b.L
                @test spmatrix(SigmaMinusPlus(i,j),b, Int) == adjoint(spmatrix(SigmaMinusPlus(j,i),b, Int))
            end
        end
        
        @testset "NumberConservedBasis($L, $N)" for (L,N) in [(3,1), (4,2), (6,4)]
            b = NumberConservedBasis(L,N)
            @testset "i=$i, j=$j" for i in 1:b.L, j in 1:b.L
                @test spmatrix(SigmaMinusPlus(i,j),b, Int) == adjoint(spmatrix(SigmaMinusPlus(j,i),b, Int))
            end
        end
    end
    
end