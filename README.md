# Optimization Algorithms

This repository includes:
1.  Nonlinear-equation solvers inside the `nonlinSolvers.jl` module. The module include the following method
    1. Bisection method: `bisection`.
    2. Fixed point method. `fixedPoint`.
    3. Newton method (for nonlinear equations). `nonlinNewton`.
    4. Secant method. `secant`.

    5. Gradient check. `checkDerivative`. This function serves as a quick check for a function and its gradient (returns `true` if the gradient function matches the numerically-estimated gradient using the given function at a random point). Similarly, it can check for the Hessian (optional).

2. Line search methods to find stationary points of functions.
    1. Global line search method. `globalLineSearch`. Takes a the functions `getSearchDirection` and `getStepSize` as arguments which allows it to be used for multiple line-search schemes such as the gradient method, Newton's method, and so on. Similiarly, the `getStepSize` function allows multiple line-search rules to be implemented such as the Armijo rule or the Wolfe-Powell rule.

    2. Gradient method. `gradientMethod`. Special case of `globalLineSearch` where `getSearchDirection = (∇f, x) ↦ -∇f(x)`.

    3. Step size functions:
        1.  Armijo rule.


The scripts are written in Julia programming language.

