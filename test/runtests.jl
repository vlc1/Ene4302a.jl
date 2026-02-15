using Plots
using Ene4302a
using Test

@testset "Sinusoidal solution" begin
    x, y = 0., ones(1)

    eq = Sinusoidal()
 #   z = eq(x, y)

    sol = Solution(x, y, eq)

    y′ = sol(x)

    step = 0.1
    fwd = ForwardEuler(step, eq)
    bwd = BackwardEuler(step, similar(y), eq)
    mid = Midpoint(step, similar(y), eq)

    p = plot(sol, 0.:0.01:1., 1)
    scatter!(p, fwd, (0, 1), x, ones(1), 1)
    scatter!(p, bwd, (0, 1), x, ones(1), 1)
    scatter!(p, mid, (0, 1), x, ones(1), 1)

    @test isa(p, Plots.Plot)
end
