using Test
include("..\\nonlinSolvers.jl");

@testset "bisection method" begin
    # without exporting
    @time @test bisection(x->x^2-1,0,5) ≈ 1 atol=0.1;
    @time @test bisection(x->(x-1)*(x+3),0,5,ϵ=1e-15) ≈ 1 atol=1e-15;
    @time @test bisection(x->(x-1)*(x+3),-5,5, ϵ=1e-15) ≈ -3 atol=1e-15;
    @time @test bisection(x->x^2-1,1,1) == 1;
    
    # export data
    @time @test bisection(x->(x-1)*(x+3),0,5,ϵ=1e-15,exportData=true) ≈ 1 atol=1e-15;
    @time @test bisection(x->(x-1)*(x+3),0,5,ϵ=1e-15,exportData=true,fileName="test1_") ≈ 1 atol=1e-15;
    @time @test bisection(x->(x-1)*(x+3),0,5,ϵ=1e-15,exportData=true,fileName="test2_",fileDir = "C:\\Users\\albaa\\Google Drive\\McGill OneDrive\\McGill\\Courses\\13_WINTER 2019\\MECH 579\\Code2\\Julia\\data\\") ≈ 1 atol=1e-15;
    try
        # test new directory
        # create new directory if not there already
        mkdir("testDir\\")  
    catch
        nothing
    end
    @time @test bisection(x->(x-1)*(x+3),0,5,ϵ=1e-15,exportData=true,fileName="test2_",fileDir = "testDir\\") ≈ 1 atol=1e-15;


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