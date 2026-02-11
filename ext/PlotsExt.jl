module PlotsExt

using Ene4302a
using Plots

Plots.@recipe function f(sol::Solution, x::AbstractVector, i)
    x, sol.(x, i)
end

Plots.@recipe function f(num::ForwardEuler, (a, b)::Tuple, i)
    (; x, y) = num

    xs = similar(y, typeof(x), 0)
    yis = similar(y, 0)

    for (x, y) in num
        x ≥ a && (push!(xs, x); push!(yis, y[i]))

        x > b && break
    end

    xs, yis
end

end
