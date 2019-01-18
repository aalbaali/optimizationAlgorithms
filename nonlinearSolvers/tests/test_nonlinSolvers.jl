using Test
include("..\\nonlinSolvers.jl");

@testset "bisection method" begin
    @time @test bisection(x->x^2-1,0,5) ≈ 1 atol=0.1;
    @time @test bisection(x->(x-1)*(x+3),0,5,ϵ=1e-15) ≈ 1 atol=1e-15;
    @time @test bisection(x->(x-1)*(x+3),-5,5, ϵ=1e-15) ≈ -3 atol=1e-15;
    @time @test bisection(x->x^2-1,1,1) == 1;
    
    # so solution
    @time @test isnan(bisection(x->x^2+1,0,5));
    @test isnan(bisection(x->x^2-1,0,0));
    
    @test_throws MethodError bisection(x->x^2-1);       # Error: missing arguments (bounds)
    @test_throws ArgumentError bisection(x->x^2-1,1,0); # Error: lower bound > upper bound 
end;


