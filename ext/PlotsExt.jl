module PlotsExt

using Ene4302a
using Plots

using Ene4302a: TimeStepper

#=
Plots.@recipe function f(sol::Solution, x::AbstractVector, i)
    x, sol.(x, i)
end

=#

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

Plots.@recipe function f(scheme::TimeStepper{1}, xs::AbstractRange, y, i=firstindex(y))
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

#=
Plots.@recipe function f(scheme::TimeStepper{1}, (a, b)::Tuple, x, y, i)
    xs = similar(y, typeof(x), 0)
    yis = similar(y, 0)

    x ≥ a && (push!(xs, x); push!(yis, y[i]))

    while x ≤ b
        x = scheme(x, y)

        x ≥ a && (push!(xs, x); push!(yis, y[i]))
    end

    xs, yis
end
=#

end
