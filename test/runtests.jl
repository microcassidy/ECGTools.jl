using ECGTools
using Test

@testset "ECGTools.jl" begin
    @testset "frame operator" begin
        X = rand(9,100)
        ℱi = InverseFrameOperator()
        ℱ = FrameOperator()
        @test ℱi * (ℱ * X) ≈ X
        @test Float64.(collect(ECGTools.left_inverse2)) == ECGTools.left_inverse
    end
end
