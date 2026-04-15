module ECGTools
import Base:*,\

# Write your package code here.
const frame_operator = Rational[1 -1 0 0 0 0 0 0 0; 0 -1 1 0 0 0 0 0 0; -1 0 1 0 0 0 0 0 0; -1//2 1 1//2 0 0 0 0 0 0; 1 -1//2 1//2 0 0 0 0 0 0; 1//2 -1//2 1 0 0 0 0 0 0; 0 2//3 -2//3 1 0 0 0 0 0; 0 2//3 -2//3 0 1 0 0 0 0; 0 2//3 -2//3 0 0 1 0 0 0; 0 2//3 -2//3 0 0 0 1 0 0; 0 2//3 -2//3 0 0 0 0 1 0; 0 2//3 -2//3 0 0 0 0 0 1]
# const FtF_inverse = [171/289 131/289 59/289 -(48/289) -(48/289) -(48/289) -(48/289) -(48/289) -(48/289);
# 131/289 195/289 79/289 -(232/867) -(232/867) -(232/867) -(232/867) -(232/867) -(232/867);
# 59/289 79/289 115/289 24/289 24/289 24/289 24/289 24/289 24/289;
# -(48/289) -(232/867) 24/289 3209/2601 608/2601 608/2601 608/2601 608/2601 608/2601;
# -(48/289) -(232/867) 24/289 608/2601 3209/2601 608/2601 608/2601 608/2601 608/2601;
# -(48/289) -(232/867) 24/289 608/2601 608/2601 3209/2601 608/2601 608/2601 608/2601;
# -(48/289) -(232/867) 24/289 608/2601 608/2601 608/2601 3209/2601 608/2601 608/2601;
# -(48/289) -(232/867) 24/289 608/2601 608/2601 608/2601 608/2601 3209/2601 608/2601;
# -(48/289) -(232/867) 24/289 608/2601 608/2601 608/2601 608/2601 608/2601 3209/2601;]


const left_inverse = [ 40/289 -(72/289) -(112/289) 75/289 135/289 79/289 0 0 0 0 0 0;
-(64/289) -(116/289) -(52/289) 169/289 73/289 47/289 0 0 0 0 0 0;
-(20/289) 36/289 56/289 107/289 77/289 105/289 0 0 0 0 0 0;
88/867 304/867 72/289 -(124/867) 8/867 116/867 1 0 0 0 0 0;
88/867 304/867 72/289 -(124/867) 8/867 116/867 0 1 0 0 0 0;
88/867 304/867 72/289 -(124/867) 8/867 116/867 0 0 1 0 0 0;
88/867 304/867 72/289 -(124/867) 8/867 116/867 0 0 0 1 0 0;
88/867 304/867 72/289 -(124/867) 8/867 116/867 0 0 0 0 1 0;
88/867 304/867 72/289 -(124/867) 8/867 116/867 0 0 0 0 0 1]

export FrameOperator
struct ForwardOperator{Float64} <: AbstractMatrix{Float64}
    X::Matrix{Float64}
    ForwardOperator() = new{Float64}(frame_operator)
end
struct InverseOperator{Float64} <: AbstractMatrix{Float64}
    X::Matrix{Float64}
    InverseOperator() = new{Float64}(left_inverse)
end
Base.size(m::InverseOperator) = size(m.X)
Base.length(m::InverseOperator) = length(m.X)
Base.getindex(m::InverseOperator,i::Integer) = m.X[i]
Base.getindex(m::InverseOperator,I::Vararg{Int,2}) = m.X[I...]

Base.size(m::ForwardOperator) = size(m.X)
Base.length(m::ForwardOperator) = length(m.X)
Base.getindex(m::ForwardOperator,i::Integer) = m.X[i]
Base.getindex(m::ForwardOperator,I::Vararg{Int,2}) = m.X[I...]

struct FrameOperator
    F::ForwardOperator
    F_l::InverseOperator
    FrameOperator() = new(ForwardOperator(),InverseOperator())
end




export \,*
\(op::FrameOperator,X::AbstractMatrix{Float64}) = copy(op.F_l) * X
*(op::FrameOperator,X::AbstractMatrix{Float64}) = op.F * X
end
