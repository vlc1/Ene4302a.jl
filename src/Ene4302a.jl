module Ene4302a

using NLsolve

export OrdinaryDifferentialEquation,
       Propagator,
       ForwardEuler,
       BackwardEuler,
       Midpoint,
       Sinusoidal

"""

    OrdinaryDifferentialEquation{N} <: Function

Abstract base type for ordinary differential equations of order `N`.

The ODE is assumed to be explicitly defined in the form:

```math
\\y ^ {(n)} (t) = f(t, y(t), \\dot{y}(t), \\ldots, y ^ {(n-1)}(t))
```

Subtypes should implement the call signatures required by the integrators.

"""
abstract type OrdinaryDifferentialEquation{N} <: Function end

const ODE = OrdinaryDifferentialEquation

"""

    Propagator{Q<:ODE}

A callable object that takes the current state of the system and computes the next state based on the underlying ODE and the chosen step and initial condition.

"""
struct Propagator{Q<:ODE}
    eq::Q
end

function (this::Propagator)(x, y::AbstractArray, tau, i)
    (; eq) = this

    propagate(eq, x, y, tau, i)
end

include("odes/sinusoidal.jl")

"""

    Integrator{N}

Abstract base type for numerical time-stepping schemes for solving first-order initial value problems.

# Type Parameter

`N::Int`: The number of previous steps required by the scheme.

- `N=1` for single-step methods (e.g., Forward Euler)
- `N=2` for two-step methods (e.g., Adams-Bashforth), etc.

"""
abstract type Integrator{N} end

include("integrators/single.jl")

end
