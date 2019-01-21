using Dates
using Logging

# logger = SimpleLogger(stdout, 0);
# global_logger(logger);

function doNothing()
    nothing
    A = 3;
    @debug "This is a function that does nothing! Debug!"
    @info "This is a function that does nothing!" a=2 A
    @error "Not really an error"
end

function bisection(func::Function, lowerLimit::Real, upperLimit::Real; ϵ::AbstractFloat= 1e-5, ϵ_f::AbstractFloat=1e-1, exportData::Bool = false,fileName::String="", fileDir::String="")
    """ This function is designed to find roots of a nonlinear, nonconvex equation 'func': Rⁿ ↦ R (func = 0) given that:
    1) A root exists in the range x∈[lowerLimit,upperLimit].
    2) func(lowerLimit)*func(upperLimit) < 0.
    
    ϵ : |x_{n+1}-x_n| < ϵ
    ϵ_f: |f(x_{n+1})| < e_f (note that e_f shouldn't be much smaller than ϵ)
    """

    if lowerLimit > upperLimit
        throw(ArgumentError("Lower limit should be less than the upper limit"))
    end
    
    # export to a file if the user wants
    if exportData
        if fileDir == ""
            fileDir = "data\\";
        end
        
        fileHandle = open(fileDir*"bisection_"*fileName*Dates.format(Dates.now(),"yyyymmddHHMM")*".txt","w");
        # fileHandle = open(fileDir*"bisection_"*fileName*Dates.format(Dates.now(),"yyyymmddHHMM")*".txt","w");
        write(fileHandle,"func: $func\tlowerLimit:$lowerLimit\tupperLimit:$upperLimit\tϵ:$ϵ\n");
        write(fileHandle,"iterations\t(upperLimit+lowerLimit)/2\tf(x)\n");
    end
    
    # upperLimit₀ = upperLimit; # first upper limit
    # lowerLimit₀ = lowerLimit; # first lower limit
    n = 0; # number of maxIterations
    r = (upperLimit+lowerLimit)/2;
    while abs(upperLimit-lowerLimit) >= ϵ
        f_r = func(r); # store f(r) to avoid calling the function twice
        if exportData
            write(fileHandle,"$n\t$r\t$f_r\n")
        end
        
        if f_r*func(lowerLimit) > 0
            lowerLimit = r;
        else
            upperLimit = r;
        end
        
        n += 1;
        
        r = (upperLimit+lowerLimit)/2;
        
    end
    
    # checking whether the solution is close to the root or not
    f_r = func(r);
    if !isapprox(f_r, 0, atol=ϵ_f)
        r = NaN;
    end
    
    if exportData
        write(fileHandle,"$n\t$r\t$f_r\n")
        close(fileHandle);
    end
    return r;
    
end;


function fixedPoint(g::Function, x₀::Real; ϵ::AbstractFloat= 1e-5, exportData::Bool = false,fileName::String="", fileDir::String="", maxIterations::Int = convert(Int,1e6))
    """ This function is designed to find the fixed points of a continuous function g(x) (i.e., find x s.t. g(x)=x.).
    
    x₀: initial point
    ϵ (epsilon): solution tolerance (|x_actualSolution-x_founSolution|)
    exportData: user chooses whether to export data or not
    fileName: description of function (optional)
    fileDir: file driectory to store the .txt file (optional)
    
    Existence and uniqueness of a fixed-point:
    1. If g∈C0[a: b] and g([a:b])⊂ [a:b], then g has at least one fixed point in [a:b]. i.e., the equation g(x)=x has at least one solution in [a:b].
    2. Suppose, in addition, that g' exists on [a:b] such that |g'(x)|< 1, ∀ x ∈ [a:b], then there is exactly one fixed point in [a:b].
    """
    
    if exportData
        if fileDir == ""
            fileDir = "data\\";
        end
        
        if fileName!=""
            fileName*="_";
        end
        
        fileHandle = open(fileDir*"fixedPoint_"*fileName*Dates.format(Dates.now(),"yyyymmddHHMM")*".txt","w");
        write(fileHandle,"g(x): $g\tx_0:$x₀\tϵ:$ϵ\tϵ:$ϵ\tmaxIterations: $maxIterations\n");
        write(fileHandle,"iterations\tx_new\tx_old\n");
    end
    
    n = 1; # iterations
    x_old = x₀;
    
    x_new = NaN;
    
    try
        x_new = g(x_old);
    catch
        throw(ArgumentError("Function could not be evaluated"));
    end
    
    
    while abs(x_new-x_old) >= ϵ
        try
            x_new, x_old = g(x_new), x_new;
        catch
            throw(ArgumentError("Function could not be evaluated"));
            # return NaN;
        end
        
        if exportData
            write(fileHandle,"$n\t$x_new\t$x_old\n");
        end
        
        n += 1;
        
        # break if max iterations reached
        if n == maxIterations
            # return NaN;
            x_new = NaN;
            break;
        end
    end
    
    if exportData
        write(fileHandle,"$n\t$x_new\t$x_old\n");
        close(fileHandle);
    end
    return x_new;
    
end;


function nonlinNewton(f::Function,∇f::Function, x₀::Union{Array,Real}; ϵ::AbstractFloat= 1e-5, maxIterations::Int = convert(Int,1e6), exportData::Bool = false,fileName::String="", fileDir::String="")
    """ Newton's method to sove nonlinear equations.
    Update scheme:
    x_{k+1} = x_k - inv(∇f(x_k))*f(x_k).
    Instead of solving for inv(∇f(x_k)), the following will be used
                            x_{k+1} = x_k + d_k.
                where d_k is obtained by solving:
                ∇f(x_k)d_{k} = -F(x_k).
                """
                
                # quick gradient check
                if !checkDerivative(f,∇f,x=x₀)
                    throw(ArgumentError("Gradient function: ∇f(x) did not conform with the given function f(x)."))
                end
                
                
                # exporting data
    if exportData
        if fileDir == ""
            fileDir = "data\\";
        end
        
        if fileName!=""
            fileName*="_";
        end
        
        fileHandle = open(fileDir*"nonlinNewton_"*fileName*Dates.format(Dates.now(),"yyyymmddHHMM")*".txt","w");
        write(fileHandle,"f(x): $f\t∇f(x): $∇f\tx_0:$x₀\tϵ:$ϵ\tϵ:$ϵ\tmaxIterations: $maxIterations\n");
        write(fileHandle,"iterations\tx_new\tf(x_new)\t∇f(x_new)\td_new\n");
    end
    
    n = 1; # iterations
    x = x₀;
    f_val = f(x);
    ∇f_val = ∇f(x);
    d = NaN;
    while abs(f_val) >= ϵ
        d = -∇f_val'\f_val; # solving for d
        x += d;
        n += 1;
        f_val = f(x);
        ∇f_val = ∇f(x);
        
        if exportData
            write(fileHandle,"$n\t$x\t$f_val\t$∇f_val\t$d\n");
        end
        
        # break if max iterations reached
        if n == maxIterations
            x = NaN;
            break;
        end
        
    end
    
    if exportData
        write(fileHandle,"$n\t$x\t$f_val\t$∇f_val\t$d\n");
        close(fileHandle);
    end
    return x;
end


function secant(f::Function, x₀::Real=rand()*-10, x₁::Real=rand()*10; ϵ::AbstractFloat= 1e-5, maxIterations::Int = convert(Int,1e6), exportData::Bool = false,fileName::String="", fileDir::String="")
    """
    The secant method is designed to solve for the roots of a nonlinear equation f(x) = 0.
    This method requires 2 intial conditions: x₀ and x₁.

    """
    if exportData
        if fileDir == ""
            fileDir = "data\\";
        end
        
        if fileName!=""
            fileName*="_";
        end
        
        fileHandle = open(fileDir*"secant_"*fileName*Dates.format(Dates.now(),"yyyymmddHHMM")*".txt","w");
        write(fileHandle,"f(x): $f\tr_l: $x₀\tr_r:$x₁\tϵ:$ϵ\tϵ:$ϵ\tmaxIterations: $maxIterations\n");
        write(fileHandle,"iterations\tx_new\tf(x)\n");
    end
    
    n = 1; # iterations
    f₀ = f(x₀); # func valued at x₀
    f₁ = f(x₁); # func valued at x₁
    
    updateIteration(f₀,f₁,r₀,r₁) = r₁ - f₁*(r₁-r₀)/(f₁-f₀);
    x = updateIteration(f₀,f₁,x₀,x₁);   # x is the latest iteration
    # while abs(x-x₀) >= ϵ
    while abs(f₁-f₀) >= ϵ
        
        if exportData
            write(fileHandle,"$n\t$x₁\t$f₁\n");
        end
        
        x₀ = x₁;
        x₁ = x;
        f₀ = f(x₀);
        f₁ = f(x₁);
        x = updateIteration(f₀,f₁,x₀,x₁);   # x is the latest iteration
        
        # break if max iterations reached
        if n == maxIterations
            x = NaN;
            break;
        end

        n+=1;
        
    end
    
    if exportData
        write(fileHandle,"$n\t$x₁\t$f₁\n");
        close(fileHandle);
    end

    if abs(f(x)) <= ϵ
        return x;
    else
        return NaN;
    end

end

function checkDerivative(f::Function,∇f::Function; ∇²f::Union{Function,Nothing}=nothing, x::Union{Real,Array},t::AbstractFloat=1e-5)
    """ f: Rⁿ↦{true,false}. This function serves as a 'gradient check'. Give it the function and gradient and it will assess it.
    """

    if typeof(x) <: Array
        d = rand(size(x)[1]); # arbitrary direction of same size as x (given that x is a column array)
    else
        d = 1;
    end
    
    f_prime = (f(x+t*d)-f(x))/t;

    if ∇²f == nothing
        ans = isapprox(f_prime, ∇f(x)'*d, atol=0.01);
        if !ans
            @warn "∇f doesn't seem to be right";
        end

        return ans;
    else
        return isapprox(f_prime, ∇f(x)'*d, atol=t*10) && checkDerivative(∇f,∇²f,x=x);
    end
end
