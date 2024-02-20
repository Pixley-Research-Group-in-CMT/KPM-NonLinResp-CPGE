using CPGE
using Test

@testset "CPGE.jl" begin
    # Write your tests here.
    L = 30
    Ls = [L,L,L]
    M = 2.0 # 2.0
    t = 1.0
    γ = 0.8  
    W = 0#sqrt(0.01)
    twists = 2*pi*rand(3) #[pi,pi,pi]
    boundary = exp.(1im.*twists)
    @time H0 = TwoBandModel(Ls, t, M, γ;boundary=boundary)
    println("Pass H0!")
end
