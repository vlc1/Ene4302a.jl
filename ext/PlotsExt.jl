module PlotsExt

using Ene4302a
using Plots

using Ene4302a: Integrator

Plots.@recipe function f(flow::Propagator, xs::AbstractRange, y, i=firstindex(y))
    n, tau = length(xs), step(xs)

    yis = similar(y, n)

    x = first(xs)
    j = firstindex(yis)
    yis[j] = y[i]

    while j < lastindex(yis)
        j = nextind(yis, j)
        x, yis[j] = flow(x, y, tau, i)
    end

    xs, yis
end

Plots.@recipe function f(scheme::Integrator{1}, xs::AbstractRange, y, i=firstindex(y))
    n, tau = length(xs), step(xs)

    yis = similar(y, n)

    x = first(xs)
    j = firstindex(yis)
    yis[j] = y[i]

    while j < lastindex(yis)
        j = nextind(yis, j)
        x = scheme(x, y, tau)
        yis[j] = y[i]
    end

    xs, yis
end

end
