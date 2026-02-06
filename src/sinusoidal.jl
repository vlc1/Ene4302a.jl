"""

 .   solution(t; y₀ = one(t), λ = -0.2, ω = 4., A = 0.5)

Analytical solution to the initial value problem with sinusoidal forcing.

# Arguments

- `t`: time
- `y₀`: initial condition
- `λ`: decay rate parameter
- `ω`: angular frequency of the forcing term
- `A`: amplitude of the forcing term

# Mathematical formulation

This function computes the analytical solution:

```math
y(t) = \\exp(-\\lambda t) y_0 + \\frac{A}{\\lambda^2 + \\omega^2} \\left( \\lambda \\sin(\\omega t) - \\omega \\cos(\\omega t) + \\omega \\exp(-\\lambda t) \\right)
```

# Returns

- The value of ``y(t)`` at time `t`

"""

solution(t; y₀ = one(t), λ = -0.2, ω = 4., A = 0.5) =
    exp(-λ * t) * y₀ + A * (λ * sin(ω * t) - ω * cos(ω * t) + ω * exp(-λ * t)) / (λ ^ 2 + ω ^ 2)

"""

 .   model(t; y₀ = one(t), λ = -0.2, ω = 4., A = 0.5)

Right-hand side of the differential equation for a sinusoidally forced first-order ODE.

# Arguments

- `t`: time
- `y₀`: initial condition
- `λ`: decay rate parameter
- `ω`: angular frequency of the forcing term
- `A`: amplitude of the forcing term

# Mathematical formulation

This function corresponds to the right-hand side of the initial value problem:

```math
\\dot{y}(t) = A \\sin(\\omega t) - \\lambda y(t)
```

with initial condition ``y_0``.

# Returns

- The time derivative ``\\dot{y}(t)`` at time `t`

"""
model(t; y₀ = one(t), λ = -0.2, ω = 4., A = 0.5) =
    A * sin(ω * t) - λ * solution(t; y₀, λ, ω, A)