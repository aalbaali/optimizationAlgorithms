# Optimization Algorithms

This repository includes:
1.  Nonlinear-equation solvers inside the `nonlinSolvers.jl` module. The module include the following method
    1. Bisection method: `bisection`.
    2. Fixed point method. `fixedPoint`.
    3. Newton method (for nonlinear equations). `nonlinNewton`.
    4. Secant method. `secant`.
    5. Gradient check. `checkDerivative`. This function serves as a quick check for a function and its gradient (returns `true` if the gradient function matches the numerically-estimated gradient using the given function at a random point). Similarly, it can check for the Hessian (optional).


The scripts are written in Julia programming language.

