# InitialValueProblems

[![Build Status](https://github.com/vlc1/InitialValueProblems.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vlc1/InitialValueProblems.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://vlc1.github.io/InitialValueProblems.jl/dev/)

Light-weight initial value problems for ordinary differential equations, with single-step integrators and a small model catalog.

## Quick Example

```julia
using Plots
using InitialValueProblems

# Define a sinusoidal ODE: dy/dx = A*sin(ω*x) - λ*y
eq = Sinusoidal(lambda = 0.2, omega = 4.0, amp = [0.5])

# Initial condition and time range
y0 = [1.0]
xs = range(0.0, 10.0, length=100)

# Plot analytical solution
prop = Propagator(eq)
plot(prop, xs, y0, 1, label="Analytical", linewidth=2)

# Compare with numerical integrators
tau = step(xs)
plot!(ForwardEuler(eq), xs, copy(y0), 1, (), label="Forward Euler")
plot!(BackwardEuler(eq), xs, copy(y0), 1, ([0.0],), label="Backward Euler")
plot!(Midpoint(eq), xs, copy(y0), 1, ([0.0],), label="Midpoint")

xlabel!("x")
ylabel!("y")
```

The package provides Plots.jl recipes for both analytical propagators and numerical integrators, making it easy to visualize and compare solutions.

