module PlotsExt

using Ene4302a
using Plots

Plots.@recipe function f(sol::Solution, x::AbstractVector, i)
    x, sol.(x, i)
end

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

end
