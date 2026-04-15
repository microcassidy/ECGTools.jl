using ECGTools
using Test

@testset "ECGTools.jl" begin
    @testset "frame operator" begin
        X = rand(9,100)
        ℱ = FrameOperator()
        @test ℱ \ (ℱ * X) ≈ X
    end
end
