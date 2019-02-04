using Test
using Logging
logger = SimpleLogger(stdout, Logging.Debug)
# logger = SimpleLogger(stdout, Logging.Warn)
# logger = SimpleLogger(stdout, Logging.Error)
global_logger(logger);

include("lineSearch.jl")


@testset "Armijo rule" begin
    f(x) = x[1]^2+(x[2]-1)^2;
    ∇f(x) = [2x[1];2(x[2]-1)];
    x = rand(2,1)
    d = -∇f(x);
    α = armijoRule(f,∇f,x,d);
    @test f(x+α*d) < f(x)
end


@testset "Wolfe_Powell rule" begin
    f(x) = x[1]^2+(x[2]-1)^2;
    ∇f(x) = [2x[1];2(x[2]-1)];
    x = rand(2,1);
    d = -∇f(x);
    α = Wolfe_Powell_rule(f,∇f,x,d, ρ=0.5);
    @test f(x+α*d) < f(x)
    
    f(x) = 2*(x.-1).^2+3;
    ∇f(x) = 4*(x-1);
    x = rand()*10;
    d = -∇f(x);
    α = Wolfe_Powell_rule(f,∇f,x,d, ρ=0.5);
    @test f(x+α*d) < f(x)
    
    # x = rand()*10;
    # d = [-∇f(x)];
    # α = Wolfe_Powell_rule(f,∇f,x,d, ρ=0.5);
    # @test f(x+α*d) < f(x);
end