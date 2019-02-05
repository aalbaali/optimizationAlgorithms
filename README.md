# Optimization Algorithms

This repository includes:
1.  Nonlinear-equation solvers inside the `nonlinSolvers.jl` module. The module include the following method
    1. Bisection method: `bisection`.
    2. Fixed point method. `fixedPoint`.
    3. Newton method (for nonlinear equations). `nonlinNewton`.
    4. Secant method. `secant`.

    5. Gradient check. `checkDerivative`. This function serves as a quick check for a function and its gradient (returns `true` if the gradient function matches the numerically-estimated gradient using the given function at a random point). Similarly, it can check for the Hessian (optional).

2. Descent methods to find stationary points of functions.
    1. Global line search method. `globalLineSearch`. Takes a the functions `getSearchDirection` and `getStepSize` as arguments which allows it to be used for multiple line-search schemes such as the gradient method, Newton's method, and so on. Similiarly, the `getStepSize` function allows multiple line-search rules to be implemented such as the Armijo rule or the Wolfe-Powell rule.

    2. Gradient method. `gradientMethod`. Special case of `globalLineSearch` where `getSearchDirection = (∇f, x) ↦ -∇f(x)`.
    3. Newton method. `exactNewton` uses the `nonlinNewton` method for solving nonlinear equations where the nonlinear system of equations to be solved is `∇²f*d = -∇f(x)`.
    4. Quasi-Newton methods. The `quasiNewton` method takes the Hessian update schemes (such as BFGS, DFP) as one of its arguments. This scheme is used in two other methods:
        1.  BFGS method. `quasiNewtonBFGS` method uses the general `quasiNewton` scheme. It uses the BFGS update.
        2.  DFP method. `quasiNewtonDFP` method uses the general `quasiNewton` scheme. It uses the DFP update.
    5. Nonlinear Conjugate Gradient methods. The `generalCG` is a general CG scheme that takes the function `updateβ` as one of its arguments. The following CG methods use this scheme.
        1. FR update. `cgFR` uses the `generalCG` scheme.
        1. PR update. `cgPR` uses the `generalCG` scheme.
        1. HS update. `cgHS` uses the `generalCG` scheme.

3. Line search methods.
    1. Armijo rule. This is a backtracking approach.
    2. Wolfe-Powell rule. This method is used in the conjugate gradient methods.


The scripts are written in Julia programming language.

