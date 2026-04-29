using ECGTools
using ECGTools:get_transform
using Test

@testset "ECGTools.jl" begin

    # "+aVL = ( I - III ) ÷2"""
    @test get_transform("aVL") == (get_transform("I") - get_transform("III"))/2
    @test -get_transform("aVR") == (get_transform("I") + get_transform("II"))/2
    @test -get_transform("aVF") == (get_transform("II") + get_transform("III"))/2
    # -aVR = ( I + II ) ÷2
    # +aVF = ( II + III) ÷2")





    @testset "frame operator" begin
        X = rand(9,100)
        ℱi = InverseFrameOperator()
        ℱ = FrameOperator()
        @test ℱi * (ℱ * X) ≈ X
        @test Float64.(collect(ECGTools.left_inverse2)) == ECGTools.left_inverse
    end
end
