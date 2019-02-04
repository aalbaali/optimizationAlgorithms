include("nonlinSolvers.jl");
include("lineSearch.jl");

using Logging;
using LinearAlgebra;
# logger = SimpleLogger(stdout, Logging.Debug);
# global_logger(logger);




# ∇f(x)
# ↦
# Rⁿ
# ∈
# dₖ
# αₖ

function generalLineSearchDescent(f::Function, ∇f::Function, x₀=1; getSearchDirection::Function, getStepSize::Union{Real,Function,Nothing}=nothing, ϵ::Real=1e-5, maxIterations::Real=1e6, exportData::Bool = false,fileName::String="", fileDir::String="")
    @debug "inside general line search. x₀:$x₀"
    """
    General line search method.
    Essentially the update scheme is:
                    x_{k+1} = x_{k} + αₖ.dₖ,
                where dₖ is the search direction and αₖ is the step size.
    The step size could be constant, or a function of the gradient (∇f: Rⁿ↦Rⁿ).
    
    getSearchDirection : ∇f,x ↦ d: A function/mapping that takes the objective function's gradient and outputs a search direction at x. Where 'f: Rⁿ ↦ R' is the objective function  and 'd ∈ Rⁿ' is the search direction.
    getStepSize        : f,∇f,x,d ↦ α. A constant or a function/mapping that takes the objective function and its gradient, and outputs a step-size:
    """
    
    @debug "before getStepSize"
    getAlpha = NaN;

    if getStepSize == nothing
        getAlpha(f,∇f,x,d) = 1; # adjust later to make it "armijo rule"
        @debug "inside getStepSize"
        
        # else if it's not a funcion , then it's a REAl number
    elseif (typeof(getStepSize) <: Real)
        @debug "getStepSize is not a function. Making a function."
        getAlpha = (f,∇f,x,d) -> getStepSize; 
        @debug "typeof(getAlpha)<:Function: $(typeof(getAlpha)<:Function)"
    else
        @debug "Keeping the same function."
        getAlpha = getStepSize;  # assign method
    end
    
    # export to a file if the user wants
    if exportData
        if fileDir == ""
            fileDir = "data\\";
        end
        
        fileHandle = open(fileDir*"generalLineSearch_"*fileName*Dates.format(Dates.now(),"yyyymmddHHMM")*".txt","w");
        write(fileHandle,"ϵ:$ϵ\n");
        write(fileHandle,"iterations\tx\tf(x)\t∇f(x)\n");
    end
    
    # quick gradient check
    xₜ = NaN; # point to test derivative at.
    if typeof(x₀) <: Array
        xₜ = rand(length(x₀),1)*10;
    else
        xₜ = rand()*10;  
    end

    if !checkDerivative(f,∇f, x=xₜ)
        @warn "∇f doesn't seem to be right";
        # nothing;
    end
    
    n = 0; # number of maxIterations
    x = x₀;
    f_val = f(x);
    ∇f_val = ∇f(x);
    d = NaN;
    α = NaN;


    @debug "Before entering the loop: " x, f_val, ∇f_val, d
    while norm(∇f_val) >= ϵ
        if exportData
            write(fileHandle,"$n\t$x\t$f_val\t$(∇f_val)\n"); # adjust
        end
        
        d = getSearchDirection(∇f,x);
        α = getAlpha(f,∇f,x,d);
        xold = x;
        x += α*d;

        xnew=x;
        n += 1;
        
        
        @debug "INSIDE loop" n, xold,f_val,∇f_val, α, d, xnew
        
        
        if n == maxIterations
            @warn "Maximum iterations reached"
            @debug "n: $n. x: $x. ∇f_val: $(∇f_val). d: $d. α: $α"
            x = NaN;
            break;
        end
        
        
        f_val = f(x);
        ∇f_val = ∇f(x);
    end
    @debug "After exiting the loop" x,f_val,∇f_val, α, d, n
    
    # checking whether the solution is close to the root or not
    if !isapprox(norm(∇f_val), 0, atol=ϵ)
        @warn "norm(∇f_val) > ϵ" n, x, ∇f_val
        x = NaN;
    end
    
    if exportData
        write(fileHandle,"$n\t$x\t$f_val\t$(∇f_val)\n"); # adjust
        close(fileHandle);
    end
    return x; # adjust 
end




# change name to steepest descent
function gradientMethod(f::Function,∇f::Function, x₀::Union{Real,Array}; ϵ::AbstractFloat= 1e-5, maxIterations::Real=1e6,getStepSize::Union{Nothing,Real,Function}=nothing, exportData::Bool = false,fileName::String="", fileDir::String="")
    """
    This function tries to find the stationary of a function using the gradient method.
    The search direction is 
    dₖ =  -∇f(xₖ).
    The current default step length is
    α = 1.
    getStepSize: f,∇f,x,d ↦ α. A constant or a function/mapping that takes the objective function and its gradient, and outputs a step-size (could be the Armijo rule or some other method)
    """
    @debug "type of x₀: $(typeof(x₀))";
    getDirection(∇f::Function, x) = -∇f(x); # gradient method (B=I)
    
    if getStepSize == nothing
        # getAlpha(f,∇f,x,d) = 1; # adjust it later to make it armijo rule
        getAlpha(f,∇f,x,d) = (f,∇f,x,d)->armijoRule(f,∇f,x,d); # let Armijo rule be the default
        # else if it's not a funcion , then it's a REAl number
    elseif (typeof(getStepSize) <: Real)
        getAlpha = (f,∇f,x,d) -> getStepSize; 
    else
        @debug "Keeping the same function."
        getAlpha = getStepSize;  # assign method/function
    end
    
    return generalLineSearchDescent(f,∇f,x₀, getSearchDirection = getDirection, getStepSize = getAlpha, ϵ = ϵ, maxIterations = maxIterations, exportData=exportData,fileName="gradient_"*fileName,fileDir=fileDir);
    
end
    
    
function exactNewton(f::Function,∇f::Function, ∇²f::Union{Matrix,Function}, x₀::Union{Real,Array}; ϵ::AbstractFloat= 1e-5, maxIterations::Int=convert(Int,1e6),getStepSize::Union{Nothing,Real,Function}=nothing, exportData::Bool = false,fileName::String="", fileDir::String="")
    """ 
    This function tries to find the stationary of a function using the Newton method.
        It simply uses the 
        The search direction is 
        dₖ =  -inv(∇²f)∇f(xₖ) 
        which is done by solving the system
        ∇²f(xₖ)*dₖ=∇f(xₖ)
        
        This method uses Newton only which may not converge! Use the globalNewton for guaranteed convergence (it uses the gradient method in conjuction with it)
        """
        
        fileName = "OptimNewton_"*fileName; 
        # use the Newton's method for solving nonlinear equations. The objective function in this case is the gradient.
        return nonlinNewton(∇f, ∇²f, x₀, ϵ=ϵ, maxIterations=maxIterations, exportData=exportData,fileName=fileName, fileDir=fileDir);
end
        
        
        
        


    function quasiNewton(f::Function,∇f::Function, x₀::Union{Real,Array};H₀= nothing, updateH::Function, ϵ::AbstractFloat= 1e-5, maxIterations::Real=1e6, exportData::Bool = false,fileName::String="", fileDir::String="")
        """ General inexact Newton method. The Hessian update function must be passed (it's passed from BFGS, or DFP) in the form updateH(H,s,y)
        """
        if typeof(x₀) != typeof(H₀*x₀)
            @warn "x₀ ($(typeof(x₀))) is not the same type as H₀*x₀ ($(typeof(H₀*x₀)))"
        end
        
        @debug "Before anything" H₀

        # H₀ must be positve definite
        if !isposdef(H₀)
            throw(ArgumentError("H₀='$H₀' is NOT positive definite"))
        end

    #initialize H₀ if not given
    if H₀ == nothing
        if typeof(x₀) <: Array
            H₀ = Matrix{Float64}(I, length(x₀), length(x₀));
        else
            H₀ = 1; # if typeof(x) <: Real, then H₀ must be real as well (it shouldn't be a matrix)
        end
    end

        # export to a file if the user wants
    if exportData
        if fileDir == ""
            fileDir = "data\\";
        end
        if fileName != ""
            fileName *= "_"
        end

        fileHandle = open(fileDir*"quasiNewton_"*fileName*Dates.format(Dates.now(),"yyyymmddHHMM")*".txt","w");
        write(fileHandle,"ϵ:$ϵ\n");
        write(fileHandle,"iterations\tx\tf(x)\t∇f(x)\tH\n");
    end

        x_old = x₀;
        fVal_old = f(x_old);
        ∇fVal_old = ∇f(x_old);
        H = H₀;
        k = 0;
        

        @debug "Before while loop" 

        while (norm(∇fVal_old)>ϵ) 
            if exportData
                write(fileHandle,"$k\t$x_old\t$fVal_old\t$(∇fVal_old)\t$H\n"); # adjust
            end
            
            @debug "Inside while loop: before anything " H ∇fVal_old
            d = -H\∇fVal_old; 
            
            @debug "Inside while loop: got direction: " d
            α = Wolfe_Powell_rule(f,∇f,x_old, d);
            # α = Wolfe_Powell_rule(f,∇f,x_old, d, α₀=1, γ=2, ρ=0.9, σ=1e-4);
            @debug "Inside while loop: got step length: " α
            
            x_new = x_old + α*d;
            s = x_new - x_old;    
            fVal_new = f(x_new);
            ∇fVal_new = ∇f(x_new);
            y = ∇fVal_new - ∇fVal_old;
            @debug "Before updating H" H, s, y
            
            H = updateH(H,s,y);
            @debug "Inside while loop: got H: " H
            
            x_old = x_new;
            ∇fVal_old = ∇fVal_new;
            k = k+1;
            
        end
        
        if exportData
            write(fileHandle,"$k\t$x_old\t$fVal_old\t$(∇fVal_old)\t$H\n"); # adjust
            close(fileHandle);
        end

        return x_old
    end
    
    function quasiNewtonBFGS()
        """ Uses the general quasiNewton scheme"""
    end
    
    function quasiNewtonDFP()
        """ Uses the general quasiNewton scheme"""
end

function generalCG()
    """ General inexact conjugate gradient method. The Hessian update function must be passed
    """
end

function cgFR()
    """Fletcher-Reeves method. Uses the general conjugate gradient scheme (generalCG)"""
end

function cgHS()
    """ Hestenes-Stiefel method. Uses the general conjugate gradient scheme (generalCG)"""
end

function cgPR()
    """ Polak-Rebiere method. Uses the general conjugate gradient scheme (generalCG)"""
end

