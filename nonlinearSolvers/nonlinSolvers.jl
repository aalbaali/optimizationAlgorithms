using Dates

function doNothing()
    println("hiii")
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
        println("Found solution: r=$r.")
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
    
    # diff = abs(g(x_new)-x_new);
    # if diff > ϵ
    #     println("g(xn)-xn: $diff");
    # end
    
    if exportData
        write(fileHandle,"$n\t$x_new\t$x_old\n");
        close(fileHandle);
    end
    return x_new;

end;



