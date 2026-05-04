using ECGTools
using ECGTools:get_transform
using Test
using LinearAlgebra

@testset "ECGTools.jl" begin

    # "+aVL = ( I - III ) ÷2"""
    @testset "operator correctness" begin
        @test get_transform("aVL") == (get_transform("I") - get_transform("III"))/2
        @test -get_transform("aVR") == (get_transform("I") + get_transform("II"))/2
        @test -get_transform("aVF") == (get_transform("II") + get_transform("III"))/2
    end


    @testset "frame operator" begin
        @test Matrix(ECGTools.F * ECGTools.Fi) |> M -> M == I
    end
end
