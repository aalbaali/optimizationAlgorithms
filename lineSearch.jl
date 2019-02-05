include("nonlinSolvers.jl");
using Logging;
using LinearAlgebra;
# logger = SimpleLogger(stdout, Logging.Debug);
# global_logger(logger);



# Armijo method 2: If ∇f is given as a function, then evaluate it at 'x' and pass the value to Amrijo method 1.
function armijoRule(f::Function,∇f::Function,x::Union{Real,Array},d::Union{Array,Real}; β::Real=0.5, σ::Real=1e-4)
    ∇f_val = ∇f(x);
    return armijoRule(f,∇f_val,x,d; β=β, σ=σ);
end

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
        # @debug "Before anything" ∇f_val,x,d,β,σ
        
        # making sure that the variables are within the bounds
        if β <= 0 || β >=1
            throw(ArgumentError("β should be in the range (0,1)"))
        elseif σ <= 0 || σ >= 1
            throw(ArgumentError("σ should be in the range (0,1)"))
        end
        
        if  ∇f_val'*d >= 0
            # throw(ArgumentError("Arguments do NOT satisfy the condition: ∇f_val'*d < 0."))
            @warn ("Arguments do NOT satisfy the condition: ∇f_val'*d < 0.")
        end
        
        
        l = 0;  # output is α:=βˡ
        f_val = f(x);
        while (f(x+(β^l)*d)>f_val+(β^l)*σ* ∇f_val'*d)
            l = l+1;
        end
        α = β^l;
    end
    
    
    # function Wolfe_Powell_rule(f::Function, ∇f_val::Vector, x::Vector, d::Vector; α₀::Real=1, γ::Real=2, ρ::Real=0.9, σ::Real=1e-4)
    #     ∇f = x-> ∇f_val;
    #     return Wolfe_Powell_rule(f,∇f,x,d, α₀=α₀, γ=γ, ρ=ρ, σ=σ);
    # end

    function Wolfe_Powell_rule(f::Function, ∇f::Function, x::Union{Real,Vector, Array}, d::Union{Real, Vector, Array}; α₀::Real=1, γ::Real=2, ρ::Real=0.9, σ::Real=1e-4, ϵ::Real=1e-5, maxIterations::Int=convert(Int,1e5))
        
        if σ <= 0 || σ >= 0.5
            throw(ArgumentError("σ '$σ' should be in the range (0,1)"))
        elseif ρ <= σ || ρ >= 1
            throw(ArgumentError("ρ '$ρ' should be in the range (σ,1)"))
        elseif α₀<=0 
            throw(ArgumentError("α₀ '$α₀' should be > 0"))
        elseif γ <=1
            throw(ArgumentError("γ '$γ' should be > 1"))
        end



    if typeof(x) != typeof(d)
        @warn "typeof(x) ($(typeof(x))) != typeof(d) ($(typeof(d)))"
    end

    α = NaN;

    if  ∇f(x)'*d >= -ϵ
        # throw(ArgumentError("Arguments do NOT satisfy the condition: ∇f_val'*d < 0."))
        @warn("Arguments do NOT satisfy the condition: ∇f_val'*d < 0.")
    end

    ϕ(α) = f(x+α*d);
    ϕ_dot(α) = ∇f(x+α*d)'*d;
    ϕ_dot₀ = ϕ_dot(0);
    ψ(α) = ϕ(α)-ϕ(0)- σ*α*ϕ_dot₀; # We want it to be ≦ 0
    
    
    a = NaN;
    b = NaN; 

    αᵢ = α₀;
    i = 0;

    # @debug "Before entering the FIRST while loop" ϕ_dot₀
    

    while(true)
        
        if ( ψ(αᵢ) >=-ϵ)
            a = 0;
            b = αᵢ;        
            break;
        elseif (ϕ_dot(αᵢ)>=ρ*ϕ_dot₀)
            α = αᵢ;
            @debug "RETURN from first loop" α ϕ_dot(α) ψ(α) f(x+α*d)-σ*α* 
            return α
        else
            αᵢ = γ*αᵢ;
        end
        
        
        i += 1;
    end
    aⱼ = a;
    bⱼ = b;
    
    j = 0;

    # @debug "Before entering the SECOND while loop" ϕ_dot₀ α typeof(α) x d α₀
    if (ψ(αᵢ) >= -ϵ)
        while(j<maxIterations)
            αⱼ = aⱼ+(bⱼ-aⱼ)/2;   #B1
            
            
            if (ψ(αⱼ) >= 0)             
                bⱼ = αⱼ;
                j+=1;
                # @debug αⱼ aⱼ bⱼ j ψ(αⱼ)
                if αⱼ < ϵ
                    # @debug "ERRRORRR α < $ϵ" αⱼ aⱼ bⱼ j
                    break
                end

                continue; #go to B1
            elseif (ϕ_dot(αⱼ) >= ρ*ϕ_dot₀)
                α = αⱼ;
                # @debug "Solution!" α
                return α;     #STOP 2
            else #psi(tj)<0 and phi_dot(tj)< 0
                aⱼ = αⱼ;
                j = j+1;
                # @debug ψ(αⱼ) ϕ_dot(αⱼ) aⱼ j bⱼ ∇f(x)
                continue; #go to B1
            end  
            
            # if j%100 == 0
            #     @info j;
            # end
            j = j+1;
        end
    # else
    #     # @warn "Something doesn't seem right";
    end
    
    
    #phi doet used in the psi function
    
    j = j+1;
    
    # @debug "Solution!" α

    return α;
    
end
    