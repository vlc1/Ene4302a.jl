"""

    TimeStepper{N}

Abstract base type for numerical time-stepping schemes for solving initial value problems.

# Type Parameter

- `N::Int`: The number of previous steps required by the scheme. `N=1` for single-step methods (e.g., Forward Euler), `N=2` for two-step methods (e.g., Adams-Bashforth), etc.

# Interface

Concrete subtypes must implement the call operator with `N+1` state arguments:

```julia
(::TimeStepper{N})(_, ::Vararg{AbstractArray,N}) where {N}
```

where the first argument is the current state vector (modified in-place) and the remaining `N` arguments are the previous state vectors.

"""
abstract type TimeStepper{N} end

"""

    (::TimeStepper{N})(y, args...)

Generic fallback for time-stepping schemes.

This method is called when a concrete subtype of `TimeStepper` does not implement its own call operator. It throws a `MethodError` to indicate that the scheme is not properly defined.

# Throws

- `MethodError`: Always thrown, indicating the scheme's recurrence relation is not implemented for the given arguments.

"""
(::TimeStepper{N})(_, ::Vararg{AbstractArray,N}) where {N} =
    throw(MethodError("The time-stepping scheme is not defined for the given arguments."))

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

    (scheme::ForwardEuler)(x, y)

Apply one step of the forward Euler scheme.

# Arguments

- `x::Number`: Current time
- `y::AbstractArray`: Current state (modified in-place by accumulating the time-scaled increment)

# Returns

The updated time `x + tau`.

# Details

Updates the state vector `y` according to:
```math
y \\leftarrow y + \\tau f(x, y)
```

where `f` is the ODE's right-hand side. The increment is computed in-place using the pre-allocated buffer `delta` to avoid allocations.

"""
function (this::ForwardEuler)(x, y::AbstractArray)
    (; step, eq) = this
    eq(y, x, step)
    x + step
end


struct BackwardEuler{T,A<:AbstractArray,Q<:ODE{1}} <: TimeStepper{1}
    step::T
    inc::A
    eq::Q
end

function (this::BackwardEuler)(x, y::AbstractArray)
    (; step, inc, eq) = this

    x += step

    fill!(inc, zero(eltype(inc)))
    sol = nlsolve(inc) do res, inc
        eq(res, x, y, inc, -step)
    end
    y .+= sol.zero

    x
end


struct Midpoint{T,A<:AbstractArray,Q<:ODE{1}} <: TimeStepper{1}
    step::T
    inc::A
    eq::Q
end

function (this::Midpoint)(x, y::AbstractArray)
    (; step, inc, eq) = this

    x += step / 2

    fill!(inc, zero(eltype(inc)))
    sol = nlsolve(inc) do res, inc
        eq(res, x, y, inc, -step, one(step) / 2)
    end
    y .+= sol.zero

    x += step / 2
end
