using ECGTools
using ECGTools:get_transform
using ECGTools:ECG_LEAD_IDXS
using Test
using LinearAlgebra
using StatsBase
using Plots



include("test_signal.jl")
@testset "ECGTools.jl" begin

    # "+aVL = ( I - III ) ÷2"""
    @testset "single signal function" begin single_signal();@test true end

    @testset "Electrode Covariance" begin
        cov_X = cov_electrode(ecg)
        @test cov_X ≈ cov_X'
        @test size(cov_X) == (8,8)
        @test rank(cov_X) == 8
    end
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
        terminals = ecg |> ecg_to_terminals
        @info "ecg: $(ecg[:,1])"
        @info "terminals: $(terminals[:,1])"

        ecghat = terminals |> terminals_to_ecg |> Matrix
        @info ecghat[:,1]
        @test ecghat[:,1] ≈ ecg[ECGTools.ECG_LEAD_IDXS,1]

    end
end
