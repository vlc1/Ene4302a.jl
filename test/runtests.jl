using Plots
using Ene4302a
using Test

@testset "Sinusoidal solution" begin
    x, y = 0., ones(1)

    eq = Sinusoidal()
    z = eq(x, y)

    sol = Solution(x, y, eq)

    y′ = sol(x)

    @test isapprox(y, y′)

    num = ForwardEuler(0.1, x, y, eq)

    p = plot(sol, 0:1, 1)
    plot!(p, num, (0, 1), 1)

    @test isa(p, Plots.Plot)
end
