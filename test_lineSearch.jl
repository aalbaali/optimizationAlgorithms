using Test
using Logging
# logger = SimpleLogger(stdout, Logging.Debug)
# logger = SimpleLogger(stdout, Logging.Warn)
logger = SimpleLogger(stdout, Logging.Error)
global_logger(logger);



# Armijo method 1.
function armijoRule(f::Function,∇f_val::Union{Real,Array},x::Union{Real,Array},d::Union{Array,Real}; β::Real=0.5, σ::Real=1e-4)
    """
    The Armijo rule is a function that gives a step length α for a given x value for a given d (given that  ∇f(x)ᵀd < 0). 
        It's a backtracking approach.
        
        α = max_{l∈N₀} {βˡ s.t. f(x+βˡd) <= f(x) + βˡσ∇f(x)ᵀd }.
        Which is equivalent to finding the minimum l (integer)
        f: Rⁿ↦ R. (Function)
        ∇f_val ∈ Rⁿ (gradient evaluated at x).
        x ∈ Rⁿ
        d ∈ Rⁿ
        β,σ ∈ (0,1) ⊂ R
        """
        @debug "Before anything" ∇f_val,x,d,β,σ
        
        # making sure that the variables are within the bounds
        if β <= 0 || β >=1
            throw(ArgumentError("β should be in the range (0,1)"))
        elseif σ <= 0 || σ >= 1
            throw(ArgumentError("σ should be in the range (0,1)"))
        end
        
        if  ∇f_val'*d >= 0
            throw(ArgumentError("Arguments do NOT satisfy the condition: ∇f_val'*d >= 0."))
        end
        
        
        l = 0;  # output is α:=βˡ
        f_val = f(x);
        while (f(x+(β^l)*d)>f_val+(β^l)*σ* ∇f_val'*d)
            l = l+1;
        end
        α = β^l;
    end
    
    # Armijo method 2: If ∇f is given as a function, then evaluate it at 'x' and pass the value to Amrijo method 1.
    function armijoRule(f::Function,∇f::Function,x::Union{Real,Array},d::Union{Array,Real}; β::Real=0.5, σ::Real=1e-4)
        ∇f_val = ∇f(x);
        return armijoRule(f,∇f_val,x,d; β=β, σ=σ);
    end