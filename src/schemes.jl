# represents numerical solution using Forward Euler method
# it would be better to think of those as "iterators" rather `<: Function`
struct ForwardEuler{T,Q<:ODE{1},C<:IC{T}} <: Function
    eq::Q
    ic::C
    tau::T
end

const FE = ForwardEuler
