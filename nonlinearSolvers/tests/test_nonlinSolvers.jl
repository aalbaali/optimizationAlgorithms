using Test
include("..\\nonlinSolvers.jl");

#=
@testset "bisection method" begin
    # without exporting
    @time @test bisection(x->x^2-1,0,5) ≈ 1 atol=0.1;
    @test bisection(x->1000x^2-1,0,5) ≈ sqrt(1/1000) atol=0.1;
    @time @test bisection(x->(x-1)*(x+3),0,5,ϵ=1e-15) ≈ 1 atol=1e-15;
    @time @test bisection(x->(x-1)*(x+3),-5,5, ϵ=1e-15) ≈ -3 atol=1e-15;
    @time @test bisection(x->x^2-1,1,1) == 1;
    
    # export data
    @time @test bisection(x->(x-1)*(x+3),0,5,ϵ=1e-15,exportData=true) ≈ 1 atol=1e-15;
    @time @test bisection(x->(x-1)*(x+3),0,5,ϵ=1e-15,exportData=true,fileName="test1") ≈ 1 atol=1e-15;
    @time @test bisection(x->(x-1)*(x+3),0,5,ϵ=1e-15,exportData=true,fileName="test2",fileDir = "C:\\Users\\albaa\\Google Drive\\McGill OneDrive\\McGill\\Courses\\13_WINTER 2019\\MECH 579\\Code2\\Julia\\data\\") ≈ 1 atol=1e-15;
    try
        # test new directory
        # create new directory if not there already
        mkdir("testDir\\")  
    catch
        nothing
    end
    @time @test bisection(x->(x-1)*(x+3),0,5,ϵ=1e-15,exportData=true,fileName="test2",fileDir = "testDir\\") ≈ 1 atol=1e-15;


    # no-solution tests
    @time @test isnan(bisection(x->x^2+1,0,5));
    @test isnan(bisection(x->x^2-1,0,0));
    
    # exceptions
    @test_throws MethodError bisection(x->x^2-1);       # Error: missing arguments (bounds)
    @test_throws ArgumentError bisection(x->x^2-1,1,0); # Error: lower bound > upper bound 
    
    # try deleting folder
    try
        rm("testDir\\",recursive=true)  
    catch
        nothing
    end

end;
=#
@testset "fixedPoint method" begin
    # without exporting
    f(x) = x^2+2x+1; # ans: -1
    g₁(x) = -(1+2x)/x; # fixed point of f
    g₂(x) = sqrt(-(1+2x));
    @time @test fixedPoint(g₁, 1) ≈ -1 atol=0.1;
    
    xSol = fixedPoint(g₁,1,ϵ=1e-15);
    @time @test f(xSol) ≈ 0 atol=1e-15;    #f(xSol) ≈ 0
    
    # plotting
    @time @test fixedPoint(g₁,1,exportData=true) ≈ -1 atol=0.1;
    @time @test fixedPoint(g₁,1,exportData=true,fileName="test",fileDir = "testDir\\") ≈ -1 atol=0.1;
    
    # no solutions
    @time @test isnan(fixedPoint(x->1/x,5,maxIterations=1000));
    @time @test isnan(fixedPoint(x->1/x,5,maxIterations=convert(Int,1e6),exportData=true));
    
    
    # Exeptions
    @test_throws ArgumentError fixedPoint(g₂, 1);
    
end
