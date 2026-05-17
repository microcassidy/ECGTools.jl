module ECGTools
using StatsBase
using LinearAlgebra
using LinearMaps

export cov_electrode

const BLOCK_SIZE = 1024
include("wfdbhelpers.jl")

# Write your package code here.
const TRANSFORM_LEAD_ORDERING = String["I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"]


const frame_operator_X = Rational{Int64}[6 -6 0 0 0 0 0 0 0; 0 -6 6 0 0 0 0 0 0; -6 0 6 0 0 0 0 0 0; -3 6 -3 0 0 0 0 0 0; 6 -3 -3 0 0 0 0 0 0; 3 3 -6 0 0 0 0 0 0; 0 4 -4 6 0 0 0 0 0; 0 4 -4 0 6 0 0 0 0; 0 4 -4 0 0 6 0 0 0; 0 4 -4 0 0 0 6 0 0; 0 4 -4 0 0 0 0 6 0; 0 4 -4 0 0 0 0 0 6]
const frame_operator_lambda = Rational{Int64}(1//6)
const ECG_LEAD_IDXS = Int64[2:3;7:12]
const ForwardOperator = frame_operator_lambda * frame_operator_X[ECG_LEAD_IDXS,:]
const ForwardOperatorPinv = pinv(ForwardOperator)
const null_space = [ones(3);zeros(6)]
const SamplingOperator = I(9)[[1;3:end],:]


export select_leads
select_leads(X::AbstractMatrix{T}) where {T} = size(X,1) == 12 ? X[ECG_LEAD_IDXS,:] : error("expected 12 leads received $(size(X,1))")


const CovarianceExchangeOperator = SamplingOperator * ForwardOperatorPinv





const frame_operator = LinearMap(frame_operator_X) * frame_operator_lambda



include("covariance_estimator.jl")






# @info nullspace
@info [frame_operator_X[ECG_LEAD_IDXS,:]*frame_operator_lambda;null_space'] |> x-> (@info "condition $(cond(x))")
@info [frame_operator_X[ECG_LEAD_IDXS,:]*frame_operator_lambda;null_space'] |> rank
T = frame_operator_X[ECG_LEAD_IDXS,:]*frame_operator_lambda
@info nullspace(T)
const forwardQR = qr([frame_operator_X[ECG_LEAD_IDXS,:]*frame_operator_lambda;null_space'])


forward_operator = frame_operator_lambda*frame_operator_X


const Ψ_terminals= Rational{Int64}[I(8) zeros(8)] #terminals







const F = frame_operator_lambda * LinearMap(Rational{Int64}[ECGTools.frame_operator_X[ECG_LEAD_IDXS,:];inv(frame_operator_lambda) * null_space'])
const Fi = Rational{Int64}(1//3) * LinearMap(Rational{Int64}[5 3 1 1 1 1 1 1 -1; 2 3 1 1 1 1 1 1 -1; 5 6 1 1 1 1 1 1 -1; 2 2 3 0 0 0 0 0 0; 2 2 0 3 0 0 0 0 0; 2 2 0 0 3 0 0 0 0; 2 2 0 0 0 3 0 0 0; 2 2 0 0 0 0 3 0 0; 2 2 0 0 0 0 0 3 0])


get_transform(s::String) = findfirst(s .== TRANSFORM_LEAD_ORDERING) |> idx -> frame_operator_X[idx,:]* frame_operator_lambda


const op_null_fill = [Matrix{Rational{Int64}}(I(9)[1:8,:]);-null_space'] |> LinearMap
const op_Ψ_terminals = LinearMap(Ψ_terminals)
const OP_FORWARD = F * op_null_fill * op_Ψ_terminals' #TERMINAL -> LEADS
const OP_INVERSE = op_Ψ_terminals * Fi #LEADS -> TERMINALS
const OP_ECG_SAMPLING = I(12)[ECG_LEAD_IDXS,:] |> Matrix{Rational{Int64}} |> LinearMap
const OP_NO_NULL = LinearMap(Matrix{Rational{Int64}}([I(8) zeros(8)]))'


#D
const M_augmented_recover = Int64[0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0; -1 -1 0 0 0 0 0 0 0 0 0 0; 2 -1 0 0 0 0 0 0 0 0 0 0; -1 2 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0]
const M_augmented_recover_scale = 1//2

"""
recover the augmented leads from lead-space signal. Operator is 12 x 12 and overwrites
rows 4:6 of the input Matrix
"""
function augmented_recover!(X::Int64)::Float64
    X += M_augmented_recover_scale*(M_augmented_recover*X)
    nothing
end








export ecg_to_terminals
function ecg_to_terminals(X::AbstractMatrix{T} ) where {T}
    # QRx = y ⟺ Rx = Q'y ⟺ x = R⁻¹Q'y
    size(X,1) == 12 || error("expected 12 channels")
    t =(UpperTriangular(forwardQR.R) \ forwardQR.Q'*[X[ECG_LEAD_IDXS,:];zeros(size(X,2))'])
    t[[1;3:end],:]
    # [1:end-1,:]
    # OP_INVERSE * OP_NO_NULL * OP_ECG_SAMPLING * LinearMap(X)
    # OP_NO_NULL * OP_ECG_SAMPLING * LinearMap(X)
    # OP_INVERSE *
end


export terminals_to_ecg
# terminals_to_ecg(X::LinearMap ) =  OP_NO_NULL' * OP_FORWARD * X
function terminals_to_ecg(X::AbstractMatrix{T}) where {T}
    # xn = nullspace
    xn = [-1 -1] * X[1:2,:]
    # tt = (forwardQR.Q * UpperTriangular(forwardQR.R) * [xn;n])[1:end-1,:]
    # I(9)[[1;3:end],:]'X
    (forwardQR.Q * UpperTriangular(forwardQR.R) * [X[1:1, begin:end];xn;X[2:end,begin:end]])[1:end-1,:]

end
# function terminals_to_ecg(X::ScaledMap{T} ) where {T}
#     G = LinearMap(Matrix{Rational}(6 * ECGTools.frame_operator)) * (1 // 6)
#     #Add the redunant terminal back in and then expand to 12-lead ecg
#     G *  * LinearMap(Ψ_terminals)'* X
# end

end
