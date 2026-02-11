"""

    ForwardEuler(tau, x, y, eq)

A stateful iterator for time-stepping an ODE using the forward Euler scheme.

Each iteration yields a tuple `(x, z)` where:

- `x`: The current time
- `z`: The current numerical solution (pre-allocated buffer)

# Fields

- `tau`: Time step size
- `x`: Current time
- `y`: Initial condition
- `z`: Solution buffer (reused across iterations)
- `dot`: Buffer for in-place model evaluation
- `eq`: The differential equation to solve

!!! warning

    This iterator reuses the same mutable buffers (`z` and `dot`) for efficiency. All yielded tuples reference the same `z` array, so its values change at each iteration. To collect the solution history, create independent copies:
    ```julia
    solver = ForwardEuler(tau, x, y, eq)
    history = [copy(z) for (x, z) in solver]
    ```

    For a specific component `i` up to time `t_max`:
    ```julia
    component_i = [z[i] for (x, z) in solver if x ≤ t_max]
    ```

"""
struct ForwardEuler{T,N,A<:AbstractArray{T,N},S,Q<:ODE{1}}
    tau::S
    x::S
    y::A
    z::A
    dot::A
    eq::Q
end

ForwardEuler(tau, x, y, eq) =
    ForwardEuler(tau, x, y, similar(y), similar(y), eq)

Base.IteratorSize(::Type{<:ForwardEuler}) = Base.IsInfinite()

Base.IteratorEltype(::Type{<:ForwardEuler}) = Base.HasEltype()

Base.eltype(::Type{<:ForwardEuler{T,N,A,S}}) where {T,N,A,S} = Tuple{S,A}

function Base.iterate(iter::ForwardEuler)
    (; x, y, z) = iter
    copy!(z, y)
    (x, z), x
end

function Base.iterate(iter::ForwardEuler, x)
    (; tau, z, dot, eq) = iter

    z .+= tau .* eq(dot, x, z)
    x += tau

    (x, z), x
end
