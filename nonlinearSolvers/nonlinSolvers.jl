

function doNothing()
    println("hiii")
end

function bisection(func::Function, lowerLimit::Real, upperLimit::Real; ϵ::AbstractFloat= 1e-5)
    """ This function is designed to find roots of a nonlinear, nonconvex equation 'func': Rⁿ ↦ R (func = 0) given that:
        1) A root exists in the range x∈[lowerLimit,upperLimit].
        2) func(lowerLimit)*func(upperLimit) < 0."""

    if lowerLimit > upperLimit
        throw(ArgumentError("Lower limit should be less than the upper limit"))
    end

    upperLimit₀ = upperLimit; # first upper limit
    lowerLimit₀ = lowerLimit; # first lower limit
    n = 0; # number of maxIterations
    r = (upperLimit+lowerLimit)/2;
    while abs(upperLimit-lowerLimit) >= ϵ
        
        if func(r)*func(lowerLimit) > 0
            lowerLimit = r;
        else
            upperLimit = r;
        end
        
        n += 1;
        
        r = (upperLimit+lowerLimit)/2;
    end

    # checking whether the solution is close to the root or not
    if isapprox(func(r),0,atol=ϵ)
        return r;
    else
        return NaN;
    end

end


