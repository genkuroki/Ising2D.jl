using Ising2D
using Test

@testset "Ising2D.jl" begin
    N = 100
    S = [ising2d!(rand_ising2d(), β_ising2d, 100) for _ in 1:N]
    @test abs(sum(energy_density_ising2d.(S))/N + 1.37) < 0.015
    @test abs(sum(abs.(magnetization_ising2d.(S)))/N - 0.184) < 0.07
end
