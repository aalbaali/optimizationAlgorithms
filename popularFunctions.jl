struct objFunc
    f::Function;
    ∇f::Function;
    ∇²f::Function;
end

# rosenbrock=objFunc(x->x^2,x->2x,x->2);