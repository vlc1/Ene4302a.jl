"""

    Sinusoidal{T,A<:AbstractArray{T}} <: ODE{1}

A first-order scalar ODE with sinusoidal forcing.

This ODE represents the differential equation:
```math
\\dot{y}(t) = A \\sin(\\omega t) - \\lambda y(t)
```

The solution exhibits a combination of exponential decay (controlled by `lambda`) and periodic forcing (controlled by `omega` and `amp`). This problem is commonly used to illustrate resonance phenomena and the interaction between transient and steady-state behavior.

# Fields

- `lambda::T`: Decay rate parameter (``\\lambda``)
- `omega::T`: Angular frequency of the sinusoidal forcing term (``\\omega``)
- `amp::A`: Amplitude of the forcing term (``A``)

# Example

```julia
eq = Sinusoidal(lambda = -0.2, omega = 4.0, amp = [0.5])
```

"""
struct Sinusoidal{T,A<:AbstractArray{T}} <: ODE{1}
    lambda::T
    omega::T
    amp::A
end

Sinusoidal(; lambda = -0.2, omega = 4., amp = ones(1) / 2) = 
    Sinusoidal(lambda, omega, amp)

"""

    (::Sinusoidal)(z, x, y)

Compute the right-hand side of the differential equation for a sinusoidally forced first-order ODE in-place.

# Arguments

- `z`: Output array to store the time derivative ``\\dot{y}(t)`` (modified in-place)
- `x`: Time at which to evaluate the derivative
- `y`: Current value of the solution

# Mathematical formulation

This function corresponds to the right-hand side of the initial value problem:

```math
\\dot{y}(t) = A \\sin(\\omega t) - \\lambda y(t)
```

where ``\\lambda`` is the decay rate, ``\\omega`` is the angular frequency, and ``A`` is the forcing amplitude.

"""
#=
function (this::Sinusoidal)(z, x, y)
    λ, ω, A = this.lambda, this.omega, this.amp

    @. z = A * sin(ω * x) - λ * y
end
=#

"""

Overwrite `y` with `f(x, y) * α + y` and return `y`.

"""
function (this::Sinusoidal)(y::AbstractArray, x::Number, α::Number)
    λ, ω, A = this.lambda, this.omega, this.amp

    @. y += (A * sin(ω * x) - λ * y) * α
end

"""

Overwrite `res` with `f(x, y + inc) α + inc` and return `inc`.

"""
function (this::Sinusoidal)(res::AbstractArray, x::Number, y::AbstractArray, inc::AbstractArray, α::Number)
    λ, ω, A = this.lambda, this.omega, this.amp

    @. res .= inc + (A * sin(ω * x) - λ * (y + inc)) * α
end

"""

Overwrite `res` with `f(x, y + β * inc) α + inc` and return `inc`.

"""
function (this::Sinusoidal)(res::AbstractArray, x::Number, y::AbstractArray, inc::AbstractArray, α::Number, β::Number)
    λ, ω, A = this.lambda, this.omega, this.amp

    @. res .= inc + (A * sin(ω * x) - λ * (y + inc * β)) * α
end

#=
"""

    solution!(y, eq::Sinusoidal, y₀, Δx)

Compute the analytical solution to the initial value problem with sinusoidal forcing in-place.

# Arguments

- `y`: Output array to store the solution (modified in-place)
- `eq`: The ODE instance containing the parameters (``\\lambda``, ``\\omega``, ``A``)
- `y₀`: Initial condition value at the initial time
- `Δx`: Time difference from the initial condition, i.e., ``\\Delta x = x - x_0``

# Mathematical formulation

This function computes the analytical solution:

```math
y(x) = \\exp(-\\lambda \\Delta x) y_0 + \\frac{A}{\\lambda^2 + \\omega^2} \\left( \\lambda \\sin(\\omega \\Delta x) - \\omega \\cos(\\omega \\Delta x) + \\omega \\exp(-\\lambda \\Delta x) \\right)
```

where ``y_0`` is the initial condition, ``\\lambda`` is the decay rate, ``\\omega`` is the angular frequency, and ``A`` is the forcing amplitude.

"""
function solution!(y, eq::Sinusoidal, y₀, Δx)
    λ, ω, A = eq.lambda, eq.omega, eq.amp

    y .= exp(-λ * Δx) * y₀ + A * (λ * sin(ω * Δx) - ω * cos(ω * Δx) + ω * exp(-λ * Δx)) / (λ ^ 2 + ω ^ 2)
end
=#

function propagate(eq::Sinusoidal, x, y, step, i)
    @boundscheck checkbounds(y, i)

    λ, ω, A = eq.lambda, eq.omega, eq.amp
    @boundscheck checkbounds(A, i)

    x += step

    x, exp(-λ * x) * y[i] + A[i] * (λ * sin(ω * x) - ω * cos(ω * x) + ω * exp(-λ * x)) / (λ ^ 2 + ω ^ 2)
end

#=

"""

    solution(eq::Sinusoidal, y₀, Δx, i)

Compute the i-th component of the analytical solution to the sinusoidally forced ODE.

# Arguments

- `eq::Sinusoidal`: The ODE instance containing parameters
- `y₀`: Initial condition vector
- `Δx`: Time difference from the initial condition
- `i`: Component index to extract

# Returns

A scalar value representing the i-th component of the solution at time ``x_0 + \\Delta x``.

# Bounds checking

Component index `i` is checked against the bounds of both `y₀` and the amplitude vector `A`. Bounds checking can be disabled with `--check-bounds=no` for performance-critical code.

"""
function solution(eq::Sinusoidal, y₀, Δx, i)
    @boundscheck checkbounds(y₀, i)

    λ, ω, A = eq.lambda, eq.omega, eq.amp
    @boundscheck checkbounds(A, i)

    exp(-λ * Δx) * y₀[i] + A[i] * (λ * sin(ω * Δx) - ω * cos(ω * Δx) + ω * exp(-λ * Δx)) / (λ ^ 2 + ω ^ 2)
end

=#