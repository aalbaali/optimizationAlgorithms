include("nonlinSolvers.jl");
include("lineSearch.jl");
include("objectiveFunctions.jl");

using Logging;
using LinearAlgebra;

 
function quadraticProgramming(Q::Array,c::Array; A::Union{Array,Nothing}=nothing, b::Union{Array,Nothing}nothing, x₀=nothing, ϵ::Real=1e-5, maxIterations::Real=1e6, exportData::Bool = false,fileName::String="", fileDir::String="") 
    """
    Solving 1/2*xᵀQx + cᵀx s.t. Ax=b.
    Where Q=Qᵀ.
    Uses linear solvers of choice.

    Basically solves [Q -Aᵀ; A  0] = [-c; b];
    """
    if !issymmetric(Q)
        @warn "Matrix Q: $Q is not symmetric. It'll be converted to a symmetric matrix.";
    end
    Q = Symmetric(Q); # More efficient to store it as a symmetric matrix.
    

    
end

