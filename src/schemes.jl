abstract type TimeStepper{N} end

"""

    ForwardEuler(tau, delta, eq)

Representation of the forward Euler recurrence relation for explicit first-order ODEs.

This is a callable object that applies one step of the forward Euler scheme to update a state vector in-place. Used in conjunction with `update!` to advance the numerical solution.

# Fields

- `tau::T`: Time step size
- `delta::AbstractArray`: Pre-allocated buffer for in-place RHS evaluation (modified during each call)
- `eq::ODE{1}`: The explicit first-order ODE to solve

# Usage

```julia
y = ones(1)
eq = Sinusoidal()
tau = 0.1
scheme = ForwardEuler(tau, similar(y), eq)

x = 0.0
for _ in 1:100
    x = update!(y, scheme, x)  # Updates y in-place, returns new time
end
```

# Details

Calling the scheme via `(scheme)(y, x)` applies the recurrence relation:
```math
y^{n+1} = y^n + \\tau f(x^n, y^n)
```

All operations are performed in-place:
- The buffer `delta` is overwritten with the RHS evaluation
- The state vector `y` is updated by accumulating the time-scaled increment
- No intermediate allocations occur

The `update!` function is provided as a convenience alias that calls the scheme and returns the updated time.

"""
struct ForwardEuler{T,Q<:ODE{1}} <: TimeStepper{1}
    step::T
    eq::Q
end
#struct ForwardEuler{T,A<:AbstractArray,Q<:ODE{1}} <: Function
#    tau::T
#    delta::A
#    eq::Q
#end


"""

    (scheme::ForwardEuler)(y, x)

Apply one step of the forward Euler scheme.

# Arguments

- `y::AbstractArray`: Current state (modified in-place by accumulating the time-scaled increment)
- `x::Number`: Current time

# Returns

The updated time `x + tau`.

# Details

Updates the state vector `y` according to:
```math
y \\leftarrow y + \\tau f(x, y)
```

where `f` is the ODE's right-hand side. The increment is computed in-place using the pre-allocated buffer `delta` to avoid allocations.

"""
function (this::ForwardEuler)(y, x)
    (; step, eq) = this
    eq(y, x, step)
    x + step
end


struct BackwardEuler{T,A<:AbstractArray,Q<:ODE{1}} <: TimeStepper{1}
    step::T
    inc::A
    eq::Q
end

function (this::BackwardEuler)(y, x)
    (; step, inc, eq) = this

    x += step

    fill!(inc, zero(eltype(inc)))
    sol = nlsolve(inc) do res, inc
        eq(res, x, y, inc, -step)
    end
    y .+= sol.zero

    x
end
