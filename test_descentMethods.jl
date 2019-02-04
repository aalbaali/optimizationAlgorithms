using Test
using Logging
logger = SimpleLogger(stdout, Logging.Debug)
logger = SimpleLogger(stdout, Logging.Warn)
# logger = SimpleLogger(stdout, Logging.Error)
global_logger(logger);

using LinearAlgebra
include("descentMethods.jl")
plotResults = false;

f_global(x) = x[1]^2+(x[2]-1)^2;
∇f_global(x) = [2x[1];2(x[2]-1)];

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


@testset "quasiNewton" begin
    f(x) = (x-2)^2-1;
    ∇f(x) = 2(x-2);
    x₀ = 1;

    updateH(H,s,y) = H + (y*y')/(y'*s)-(H*s*s'*H)/(s'*H*s);
    
    # quasiNewton(f,∇f, x₀,H₀= Matrix{Float64}(I, length(x₀), length(x₀)), updateH=updateH)
    @test ∇f(quasiNewton(f,∇f, x₀, H₀= 1, updateH = updateH, exportData = plotResults, fileName = "test")) ≈ 0 atol=0.001
end


@testset "BFGS" begin
    f(x) = (x-2)^2-1;
    ∇f(x) = 2(x-2);
    x₀ = 1;

    # quasiNewton(f,∇f, x₀,H₀= Matrix{Float64}(I, length(x₀), length(x₀)), updateH=updateH)
    @test ∇f(quasiNewtonBFGS(f,∇f, x₀, H₀= 1, exportData = plotResults, fileName = "test")) ≈ 0 atol=0.001
    @test norm(∇f_global(quasiNewtonBFGS(f_global,∇f_global, [1,2], H₀= [1 0;0 1], exportData = plotResults, fileName = "test"))) ≈ 0 atol=0.001
end

@testset "DFP" begin
f(x) = (x-2)^2-1;
∇f(x) = 2(x-2);
x₀ = 1;

# quasiNewton(f,∇f, x₀,H₀= Matrix{Float64}(I, length(x₀), length(x₀)), updateH=updateH)
    @test ∇f(quasiNewtonDFP(f,∇f, x₀, H₀= 1, exportData = plotResults, fileName = "test")) ≈ 0 atol=0.001
    @test norm(∇f_global(quasiNewtonDFP(f_global,∇f_global, [1,2], H₀= [1 0;0 1], exportData = plotResults, fileName = "test"))) ≈ 0 atol=0.001
end
