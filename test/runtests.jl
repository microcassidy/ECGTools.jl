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

    @testset "lead <-> terminal conversion" begin
        f() = Int64.(rand(Int16,12,10))
        ecg = Int64[195 195 195 195 176 166 156 146 137 127; 107 107 107 107 93 88 83 78 68 68; -88 -88 -88 -88 -83 -78 -73 -68 -68 -59; -151 -151 -151 -151 -132 -127 -117 -112 -102 -98; 142 142 142 142 127 122 112 107 102 93; 10 10 10 10 5 5 5 5 0 5; -29 -29 -29 -29 -24 -20 -15 -10 0 0; 88 88 88 88 83 78 73 68 59 59; 68 68 68 68 49 39 29 20 5 -10; 293 293 293 293 278 264 249 234 210 195; 273 273 273 273 259 244 229 215 200 185; 342 342 342 342 332 322 312 303 283 273]
        terminals = ecg |> ecg_to_terminals

        ecghat = terminals |> terminals_to_ecg |> Matrix
        @test ecghat[:,1] == ecg[ECGTools.ECG_LEAD_IDXS,1]
    end
end
