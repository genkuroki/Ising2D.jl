using Ising2D
using Random: seed!
using Test

@testset "Ising2D.jl" begin
    N = 100

    seed!(4649373)
    S_ifelse = [ising2d!(IfElse(), rand_ising2d(), β_ising2d, 100) for _ in 1:N]
    @test abs(sum(energy_density_ising2d.(S_ifelse))/N + 1.37) < 0.015
    @test abs(sum(abs.(magnetization_ising2d.(S_ifelse)))/N - 0.184) < 0.07

    seed!(4649373)
    S_multifor = [ising2d!(MultiFor(), rand_ising2d(), β_ising2d, 100) for _ in 1:N]
    @test abs(sum(energy_density_ising2d.(S_multifor))/N + 1.37) < 0.015
    @test abs(sum(abs.(magnetization_ising2d.(S_multifor)))/N - 0.184) < 0.07

    @test S_ifelse == S_multifor
end
