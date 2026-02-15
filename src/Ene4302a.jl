module Ene4302a

using NLsolve

export OrdinaryDifferentialEquation,
       Propagator,
#       Solution,
#       solution,
#       solution!,
#       update!,
       ForwardEuler,
       BackwardEuler,
       Midpoint,
#       Series,
       Sinusoidal

"""

    OrdinaryDifferentialEquation{N} <: Function

Abstract base type for ordinary differential equations of order `N`.

The ODE is assumed to be explicitly defined in the form:

```math
\\y ^ {(n)} (t) = f(t, y(t), \\dot{y}(t), \\ldots, y ^ {(n-1)}(t))
```

Subtypes should implement the call signature `(::ODE)(z, x, y)` to compute the right-hand side of the ODE in-place, where `z` stores the derivative, `x` is the independent variable, and `y` is the current state.

"""
abstract type OrdinaryDifferentialEquation{N} <: Function end

const ODE = OrdinaryDifferentialEquation

(this::ODE{N})(x, args::Vararg{Any,N}) where {N} =
    (z = similar(args[1]); this(z, x, args...))

(this::ODE{N})(z, x, args::Vararg{Any,N}) where {N} =
    throw(MethodError("The right-hand side function is not defined."))

#=
"""

    solution!(y, eq::ODE, y₀, Δx)

Compute the analytical solution to an ODE in-place. Generic fallback throws an error.

Subtypes of `ODE` can implement this method to provide analytical solutions when available.

"""
solution!(_, ::ODE, _, _) = throw(ArgumentError("No analytical solution available."))

"""

    solution(eq::ODE, y₀, Δx, i)

Compute a specific component of the analytical solution to an ODE. Generic fallback throws an error.

# Arguments

- `eq::ODE`: The ODE instance
- `y₀`: Initial condition vector
- `Δx`: Time difference from the initial condition
- `i`: Component index to extract from the solution

Subtypes of `ODE` should implement this method to provide analytical solutions where available, with automatic bounds checking on the component index.

"""
solution(::ODE, _, _, _) = throw(ArgumentError("No analytical solution available."))


"""
    Solution{Q<:ODE, C<:State} <: Function

A callable object representing the analytical solution to an ODE with given initial conditions.

# Fields
- `eq::ODE`: The ordinary differential equation
- `ic::State`: The initial condition

# Usage
    sol = Solution(eq, ic)
    y = sol(x)          # Allocating: compute solution at x
    sol(y, x)           # In-place: store solution at x in y
"""
struct Solution{T,A<:AbstractArray,Q<:ODE} <: Function
    x::T
    y::A
    eq::Q
end

function (this::Solution)(x′::Number, i)
    (; x, y, eq) = this

    solution(eq, y, x′-x, i)
end

(this::Solution)(x::Number) = (y = similar(this.y); this(y, x))

function (this::Solution)(y′::AbstractArray, x′::Number)
    (; x, y, eq) = this

    solution!(y′, eq, y, x′-x)
end
=#

struct Propagator{Q<:ODE}
    eq::Q
end

function (this::Propagator)(x, y::AbstractArray, tau, i)
    (; eq) = this

    propagate(eq, x, y, tau, i)
end

include("odes.jl")
include("timestepper.jl")

end
