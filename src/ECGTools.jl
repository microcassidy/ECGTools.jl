module ECGTools
using LinearAlgebra
using LinearMaps

# Write your package code here.
const TRANSFORM_LEAD_ORDERING = String["I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"]


const frame_operator_X = Rational{Int64}[6 -6 0 0 0 0 0 0 0; 0 -6 6 0 0 0 0 0 0; -6 0 6 0 0 0 0 0 0; -3 6 -3 0 0 0 0 0 0; 6 -3 -3 0 0 0 0 0 0; 3 3 -6 0 0 0 0 0 0; 0 4 -4 6 0 0 0 0 0; 0 4 -4 0 6 0 0 0 0; 0 4 -4 0 0 6 0 0 0; 0 4 -4 0 0 0 6 0 0; 0 4 -4 0 0 0 0 6 0; 0 4 -4 0 0 0 0 0 6]
const frame_operator_lambda = Rational{Int64}(1//6)
const frame_operator = LinearMap(frame_operator_X) * frame_operator_lambda
const ECG_LEAD_IDXS = Int64[[1,3];7:12]


const Ψ_terminals= Rational{Int64}[I(8) zeros(8)] #terminals
const null_space = Rational{Int64}[-ones(3);ones(6)] #null space


const F = frame_operator_lambda * LinearMap(Rational{Int64}[ECGTools.frame_operator_X[ECG_LEAD_IDXS,:];inv(frame_operator_lambda) * null_space'])
const Fi = Rational{Int64}(1//3) * LinearMap(Rational{Int64}[5 3 1 1 1 1 1 1 -1; 2 3 1 1 1 1 1 1 -1; 5 6 1 1 1 1 1 1 -1; 2 2 3 0 0 0 0 0 0; 2 2 0 3 0 0 0 0 0; 2 2 0 0 3 0 0 0 0; 2 2 0 0 0 3 0 0 0; 2 2 0 0 0 0 3 0 0; 2 2 0 0 0 0 0 3 0])


get_transform(s::String) = findfirst(s .== TRANSFORM_LEAD_ORDERING) |> idx -> frame_operator_X[idx,:]* frame_operator_lambda
end
