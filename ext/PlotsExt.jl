module PlotsExt

using Ene4302a
using Plots

Plots.@recipe function f(sol::Solution, x::AbstractVector, i)
    x, sol.(x, i)
end

Plots.@recipe function f(sol::ForwardEuler, (a, b)::Tuple, i)
    (; eq, ic, tau) = sol
    x₀, y₀ = ic.x, ic.y

    abscissas = similar(y₀, typeof(x₀), 0)
    components = similar(y₀, 0)

    x, y, z = x₀, deepcopy(y₀), similar(y₀)
    x ≥ a && (push!(abscissas, x); push!(components, y[i]))

    while x < b
        y .+= tau .* eq(z, x, y)
        x += tau

        x ≥ a && (push!(abscissas, x); push!(components, y[i]))
    end

    abscissas, components
end

end
