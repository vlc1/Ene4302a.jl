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

Overwrite `y` with `f(x, y) * α + y` and return `y`.

"""
function (this::Sinusoidal)(x, y::AbstractArray, α::Number)
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


"""

    propagate(eq::Sinusoidal, x, y, tau, i)

Compute the `i`th component of the analytical solution to the initial value problem with sinusoidal forcing.

# Arguments

- `eq`: The ODE instance containing the parameters (``\\lambda``, ``\\omega``, ``A``)
- `x`: Initial time
- `y`: Initial condition value at the initial time
- `tau`: Time difference from the initial condition
- `i`: Index of the component to compute (must be within bounds of `y` and `A`)

# Returns

A `Tuple` containing:

- The initial time `x`
- The `i`th component of the analytical solution at time `x + tau`

# Mathematical formulation

The analytical solution is given by:

```math
y(x) = \\exp(-\\lambda \\Delta x) y_0 + \\frac{A}{\\lambda^2 + \\omega^2} \\left( \\lambda \\sin(\\omega \\Delta x) - \\omega \\cos(\\omega \\Delta x) + \\omega \\exp(-\\lambda \\Delta x) \\right)
```

where ``y_0`` is the initial condition, ``\\lambda`` is the decay rate, ``\\omega`` is the angular frequency, and ``A`` is the forcing amplitude.

"""
function propagate(eq::Sinusoidal, x, y, tau, i)
    @boundscheck checkbounds(y, i)

    λ, ω, A = eq.lambda, eq.omega, eq.amp
    @boundscheck checkbounds(A, i)

    x += tau

    x, exp(-λ * x) * y[i] + A[i] * (λ * sin(ω * x) - ω * cos(ω * x) + ω * exp(-λ * x)) / (λ ^ 2 + ω ^ 2)
end
