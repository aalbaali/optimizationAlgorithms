using Dates

function doNothing()
    println("hiii")
end

function bisection(func::Function, lowerLimit::Real, upperLimit::Real; ϵ::AbstractFloat= 1e-5, exportData::Bool = false,fileName::String="", fileDir::String="")
    """ This function is designed to find roots of a nonlinear, nonconvex equation 'func': Rⁿ ↦ R (func = 0) given that:
    1) A root exists in the range x∈[lowerLimit,upperLimit].
    2) func(lowerLimit)*func(upperLimit) < 0."""
    
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
    if !isapprox(f_r, 0, atol=ϵ)
        r = NaN;
    end
    
    if exportData
        write(fileHandle,"$n\t$r\t$f_r\n")
    end
    return r;
    
end


