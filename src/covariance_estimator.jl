"""
Y - ecg signal
"""
function cov_electrode(Y::AbstractMatrix{T}) where {T}
    size(Y,1) |> N -> N == 12 || error("expecting 12 channels in first dimension, recieved $(N)")
    C = cov(select_leads(Y),dims=2)
    Hermitian(CovarianceExchangeOperator * C * CovarianceExchangeOperator')
end
