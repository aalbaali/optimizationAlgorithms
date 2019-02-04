using Test
using Logging
# logger = SimpleLogger(stdout, Logging.Debug)
# logger = SimpleLogger(stdout, Logging.Warn)
logger = SimpleLogger(stdout, Logging.Error)
global_logger(logger);

include("descentMethods.jl")


@testset "generalLineSearchDescent" begin
    f(x) = (x-2)^2-1;
    ∇f(x) = 2(x-2);
    α = 0.5;
    getDk(∇f::Function, x) = -∇f(x);
    

    @test generalLineSearchDescent(f,∇f,2.0,getSearchDirection=getDk, getStepSize=α) ≈ 2;
    @test ∇f(generalLineSearchDescent(f,∇f,2.0,getSearchDirection=getDk, getStepSize=α, ϵ=1e-10)) ≈ 0 atol=1e-10;
    
    f(x) = x^2-1;
    ∇f(x) = 2x;
    α = 1; # keeps oscillating between two solutions
    @test isnan(generalLineSearchDescent(f,∇f,2.0,getSearchDirection=getDk, getStepSize=α));
    α = 0.5;
    @test ∇f(generalLineSearchDescent(f,∇f,2.0,getSearchDirection=getDk, getStepSize=α, ϵ=1e-10)) ≈ 0 atol=1e-10;
    ∇f(x) = 2x-1; #wrong gradient: should give a warning!
    # @test_logs Logging.Warn generalLineSearchDescent(f,∇f,2.0,getSearchDirection=getDk, getStepSize=α);
    @test_logs (:warn,"∇f doesn't seem to be right") generalLineSearchDescent(f,∇f,2.0,getSearchDirection=getDk, getStepSize=α);
    
    
    f(x) = x[1]^2+(x[2]-1)^2;
    ∇f(x) = [2x[1];2(x[2]-1)];
    @test generalLineSearchDescent(f,∇f,[2.0,2.0],getSearchDirection=getDk, getStepSize=α, ϵ=1e-10) ≈ [0,1] atol=1e-10;
    
end

@testset "gradientMethod" begin
    f(x) = (x-2)^2-1;
    ∇f(x) = 2(x-2);
    α = 0.5;
    # gradientMethod(f,\nabl)
    @test gradientMethod(f,∇f,5, getStepSize=α, ϵ=1e-10) ≈ 2 atol=1e-10;
    # gradientMethod(f,∇f,5, getStepSize=α, ϵ=1e-10) ≈ [0,1] atol=1e-10;

    
    f(x) = (x-2)^2-1;
    ∇f(x) = 2(x-2);
    α = 0.5;
    
    @test gradientMethod(f,∇f,2.0, getStepSize=α) ≈ 2;
    @test ∇f(gradientMethod(f,∇f,2.0,getStepSize=α, ϵ=1e-10)) ≈ 0 atol=1e-10;
    
    f(x) = x^2-1;
    ∇f(x) = 2x;
    α = 1; # keeps oscillating between two solutions
    @test isnan(gradientMethod(f,∇f,2.0, getStepSize=α));
    α = 0.5;
    @test ∇f(gradientMethod(f,∇f,2.0, getStepSize=α, ϵ=1e-10)) ≈ 0 atol=1e-10;
    ∇f(x) = 2x-1; #wrong gradient: should give a warning!
    # @test_logs Logging.Warn gradientMethod(f,∇f,2.0,getSearchDirection=getDk, getStepSize=α);
    @test_logs (:warn,"∇f doesn't seem to be right") gradientMethod(f,∇f,2.0,getStepSize=α);
    
    
    f(x) = x[1]^2+(x[2]-1)^2;
    ∇f(x) = [2x[1];2(x[2]-1)];
    @test gradientMethod(f,∇f,[2.0,2.0],getStepSize=α, ϵ=1e-10) ≈ [0,1] atol=1e-10;
    
end


@testset "Armijo rule" begin
    f(x) = x[1]^2+(x[2]-1)^2;
    ∇f(x) = [2x[1];2(x[2]-1)];
    @test gradientMethod(f,∇f,[2.0,2.0],getStepSize=armijoRule, ϵ=1e-10) ≈ [0,1] atol=1e-10;
end