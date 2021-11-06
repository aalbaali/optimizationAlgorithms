include("nonlinSolvers.jl");

struct functionAndItsDerivatives
    f::Function
    ∇f::Function
     ∇²f::Function
end


f(x) = (1-x[1])^2+100*(x[2]-x[1]^2)^2;
∇f(x) = [-2*(1-x[1])-400*x[1]*(x[2]-x[1]^2);200*(x[2]-x[1]^2)];
∇²f(x) = [-400*(x[2]-x[1]^2)+800*x[1]^2+2 -400*x[1]; -400*x[1] 200];
Rosenbrock = functionAndItsDerivatives(f,∇f, ∇²f)

