@testset "Pauli matrices" begin
    @test spmatrix(SigmaX(1)) == [0 1; 1 0]
    @test spmatrix(SigmaY(1)) == [0 -im; im 0]
    @test spmatrix(SigmaZ(1)) == [1 0; 0 -1]
    @test spmatrix(SigmaPlus(1)) == [0 1; 0 0]
    @test spmatrix(SigmaMinus(1)) == [0 0; 1 0]

    b = TensorBasis(1)
    @test spmatrix(SigmaX(1)) == spmatrix(SigmaX(1), b)
    @test spmatrix(SigmaY(1)) == spmatrix(SigmaY(1), b)
    @test spmatrix(SigmaZ(1)) == spmatrix(SigmaZ(1), b)
    @test spmatrix(SigmaPlus(1)) == spmatrix(SigmaPlus(1), b)
    @test spmatrix(SigmaMinus(1)) == spmatrix(SigmaMinus(1), b)

    @test spmatrix(SigmaPlus(1)) == 1/2 * (spmatrix(SigmaX(1)) + im * spmatrix(SigmaY(1)))
    @test spmatrix(SigmaMinus(1)) == 1/2 * (spmatrix(SigmaX(1)) - im * spmatrix(SigmaY(1)))

end


@testset "SigmaPlusMinus" begin

    @testset "NumberConservedBasis" begin
        
        M = spmatrix(SigmaPlusMinus(1,2), NumberConservedBasis(4,2), Int)
        M_corr = spzeros(Int, 6,6)
        M_corr[4,2] = 1
        M_corr[5,3] = 1
        @test M == M_corr
        
    end

    @testset "TensorBasis" begin
        
        M = spmatrix(SigmaPlusMinus(1,2), TensorBasis(2), Int)
        M_corr = spzeros(Int, 4,4)
        M_corr[3,2] = 1
        @test M == M_corr

        M = spmatrix(SigmaPlusMinus(2,1), TensorBasis(2), Int)
        M_corr = spzeros(Int, 4,4)
        M_corr[2,3] = 1
        @test M == M_corr

        M = spmatrix(SigmaPlusMinus(1,2), TensorBasis(3), Int)
        M_corr = spzeros(Int, 8,8)
        M_corr[5,3] = 1
        M_corr[6,4] = 1
        @test M == M_corr

        M = spmatrix(SigmaPlusMinus(2,1), TensorBasis(3), Int)
        M_corr = spzeros(Int, 8,8)
        M_corr[3,5] = 1
        M_corr[4,6] = 1
        @test M == M_corr

        M = spmatrix(SigmaPlusMinus(1,3), TensorBasis(3), Int)
        M_corr = spzeros(Int, 8,8)
        M_corr[5,2] = 1
        M_corr[7,4] = 1
        @test M == M_corr

        M = spmatrix(SigmaPlusMinus(3,1), TensorBasis(3), Int)
        M_corr = spzeros(Int, 8,8)
        M_corr[2,5] = 1
        M_corr[4,7] = 1
        @test M == M_corr

        M = spmatrix(SigmaPlusMinus(2,3), TensorBasis(3), Int)
        M_corr = spzeros(Int, 8,8)
        M_corr[3,2] = 1
        M_corr[7,6] = 1
        @test M == M_corr

        M = spmatrix(SigmaPlusMinus(3,2), TensorBasis(3), Int)
        M_corr = spzeros(Int, 8,8)
        M_corr[2,3] = 1
        M_corr[6,7] = 1
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