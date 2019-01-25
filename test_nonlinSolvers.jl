using Test
include("nonlinSolvers.jl");


@testset "checkDerivative" begin
    # check without hessian
    @test checkDerivative(x->x^2, x->2x,x=rand());
    @test checkDerivative(x->x^2, x->2x,x=rand()*10);
    @test !checkDerivative(x->x^2, x->2x+1,x=rand()*10);
    @test !checkDerivative(x->x^2, x->x+1,x=rand()*10);

    # checking with hessian
    f(x) = x[1]^2+x[1]*x[2]+x[2]^2+2;
    ∇f(x) = [2x[1]+x[2]; 2x[2]+x[1]]; 
    ∇²f(x) = [2 1;1 2];
    @test checkDerivative(f,∇f,∇²f=∇²f,x=rand(2,1));
    
    @test checkDerivative(x->x[1]^2+x[2]^2, x->[2x[1], 2x[2]],x=rand(2,1));
    @test !checkDerivative(x->x[1]^2+2x[2]^2, x->[2x[1], 2x[2]],x=rand(2,1));
    
    f(x) = x[1]^2+2x[2]^2+3x[3]*x[2];
    ∇f(x) = [2x[1], 4x[2]+3x[3], 3x[2]];
    @test checkDerivative(f,∇f,x=rand(3,1));
    @test !checkDerivative(f,x-> ∇f(x)+[0,0,x[3]],x=rand(3,1)*10);

end


@testset "bisection method" begin
    # without exporting
    @test bisection(x->x^2-1,0,5) ≈ 1 atol=0.1;
    @test bisection(x->1000x^2-1,0,5) ≈ sqrt(1/1000) atol=0.1;
    @test bisection(x->(x-1)*(x+3),0,5,ϵ=1e-15) ≈ 1 atol=1e-15;
    @test bisection(x->(x-1)*(x+3),-5,5, ϵ=1e-15) ≈ -3 atol=1e-15;
    @test bisection(x->x^2-1,1,1) == 1;
    
    # export data
    #  @test bisection(x->(x-1)*(x+3),0,5,ϵ=1e-15,exportData=true) ≈ 1 atol=1e-15;
    #  @test bisection(x->(x-1)*(x+3),0,5,ϵ=1e-15,exportData=true,fileName="test1") ≈ 1 atol=1e-15;
    # @test bisection(x->(x-1)*(x+3),0,5,ϵ=1e-15,exportData=true,fileName="test2",fileDir = "C:\\Users\\albaa\\Google Drive\\McGill OneDrive\\McGill\\Courses\\13_WINTER 2019\\MECH 579\\Code2\\Julia\\data\\") ≈ 1 atol=1e-15;
    try
        # test new directory
        # create new directory if not there already
        mkdir("testDir\\")  
    catch
        nothing
    end
    # @test bisection(x->(x-1)*(x+3),0,5,ϵ=1e-15,exportData=true,fileName="test2",fileDir = "testDir\\") ≈ 1 atol=1e-15;


    # no-solution tests
    @test isnan(bisection(x->x^2+1,0,5));
    @test isnan(bisection(x->x^2-1,0,0));
    
    # exceptions
    @test_throws MethodError bisection(x->x^2-1);       # Error: missing arguments (bounds)
    @test_throws ArgumentError bisection(x->x^2-1,1,0); # Error: lower bound > upper bound 
    
    ## try deleting folder
    # try
    #     rm("testDir\\",recursive=true)  
    # catch
    #     nothing
    # end

end;


@testset "fixedPoint method" begin
    # without exporting
    f(x) = x^2+2x+1; # ans: -1
    g₁(x) = -(1+2x)/x; # fixed point of f
    g₂(x) = sqrt(-(1+2x));
    @test fixedPoint(g₁, 1) ≈ -1 atol=0.1;
    
    xSol = fixedPoint(g₁,1,ϵ=1e-10);
    @test f(xSol) ≈ 0 atol=1e-10;    #f(xSol) ≈ 0
    
    # plotting
    #  @test fixedPoint(g₁,1,exportData=true) ≈ -1 atol=0.1;
    #  @test fixedPoint(g₁,1,exportData=true,fileName="test",fileDir = "testDir\\") ≈ -1 atol=0.1;
    
    # no solutions
    @test isnan(fixedPoint(x->1/x,5,maxIterations=1000));
    # @test isnan(fixedPoint(x->1/x,5,maxIterations=convert(Int,1e6),exportData=true));
    
    
    # Exeptions
    @test_throws ArgumentError fixedPoint(g₂, 1);
    
end

# add a test of a system of eqns
@testset "Nonlin Newton method" begin
    # without exporting
    f₁(x) = x^2+2x+1; # ans: -1
    ∇f₁(x) = 2x+2;
    ∇f₁_wrong(x) = 2x+1; # wrong gradient

    @test nonlinNewton(f₁,∇f₁,1) ≈ -1 atol=0.1;
    @test f₁(nonlinNewton(f₁,∇f₁,1,ϵ=1e-10)) ≈ 0 atol=1e-10;
    # Exeptions
    @test_throws ArgumentError nonlinNewton(f₁,∇f₁_wrong,1); # wrong gradient
    
    f₂(x) = x[1]^2 + x[2]^2 -1;
    ∇f₂(x) = [2x[1]; 2x[2]];
    @test f₂(nonlinNewton(f₂,∇f₂,[1,1],ϵ=1e-10)) ≈ 0 atol=1e-10;
    
    
    # plotting
    #  @test nonlinNewton(f₁,∇f₁,1,exportData=true) ≈ -1 atol=0.1;
    @test nonlinNewton(f₁,∇f₁,1,exportData=true,fileName="test2",fileDir = "testDir\\") ≈ -1 atol=0.1;
    
end;

@testset "secant" begin
    f(x) = x^2+2x+1;
    g(x) = x^4-16;
    @test f(secant(f)) ≈ 0 atol=0.01;
    @test f(secant(f,0,1, ϵ=1e-10)) ≈ 0 atol=1e-10;
    @test g(secant(g,0,1, ϵ=1e-10)) ≈ 0 atol=1e-10
    
    # plotting
    @test f(secant(f,exportData=true)) ≈ 0 atol=0.01;
end;
