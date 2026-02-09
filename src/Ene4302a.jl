module Ene4302a

export OrdinaryDifferentialEquation,
       InitialCondition,
       Solution,
       solution,
       solution!,
       ForwardEuler,
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

    InitialCondition{T,N,A<:AbstractArray{T,N}}

Stores the initial condition for an ODE at position/time `x` with state `y`.

# Fields

- `x::T`: Initial position or time
- `y::A`: Initial state (array with `N` dimensions)

# Constructor

    InitialCondition(x, y)
    InitialCondition(y)  # assumes x = zero(y)

"""
struct InitialCondition{T,N,A<:AbstractArray{T,N}}
    x::T
    y::A
end

const IC = InitialCondition

IC(y) = IC(zero(y), y)

"""
    Solution{Q<:ODE, C<:InitialCondition} <: Function

A callable object representing the analytical solution to an ODE with given initial conditions.

# Fields
- `eq::ODE`: The ordinary differential equation
- `ic::InitialCondition`: The initial condition

# Usage
    sol = Solution(eq, ic)
    y = sol(x)          # Allocating: compute solution at x
    sol(y, x)           # In-place: store solution at x in y
"""
struct Solution{Q<:ODE,C<:IC} <: Function
    eq::Q
    ic::C
end

function (this::Solution)(x::Number, i)
    (; eq, ic) = this
    x₀, y₀ = ic.x, ic.y

    solution(eq, y₀, x-x₀, i)
end

(this::Solution)(x::Number) = (y = similar(this.ic.y); this(y, x))

function (this::Solution)(y::AbstractArray, x::Number)
    (; eq, ic) = this
    x₀, y₀ = ic.x, ic.y

    solution!(y, eq, y₀, x-x₀)
end

include("odes.jl")
include("schemes.jl")

end
