using Plots
using InitialValueProblems
using Test

@testset "Sinusoidal solution" begin
    x, y, eq = 0., ones(1), Sinusoidal()

    ref = Propagator(eq)
    fwd = ForwardEuler(eq)
    bwd = BackwardEuler(eq, similar(y))
    mid = Midpoint(eq, similar(y))

    xs = 0.:0.1:1.

    p = plot(ref, xs, ones(1), 1)
    scatter!(p, fwd, xs, ones(1), 1)
    scatter!(p, bwd, xs, ones(1), 1)
    scatter!(p, mid, xs, ones(1), 1)

    @test isa(p, Plots.Plot)
end
