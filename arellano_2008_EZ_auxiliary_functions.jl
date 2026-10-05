util(x) = x^(1-ρ) / ((1-ρ)*(1-β))
util_inv(u) = ((1-ρ)*(1-β)*u)^(1/(1-ρ))

h = (y) -> y - max(0.0, y_const*y + y_curv*y^2)

function ex_ante_value(a::Float64, b::Float64)
    """
    Computes the log-sum-exponential of two arguments in a way that avoids numerical overflow or underflow.
    """
    vmax = max(a, b)
    return vmax + α * log(exp((a - vmax)/α) + exp((b - vmax)/α)) - α*log(2)
end

function certainty_equivalent(values, weights)
    """
    Computes the certainty equivalent in a numerically safe way that is robust
    to different EIS and risk aversion parameters
    """
    v̄ = maximum(values)
    v̄ > 0 || error("certainty_equivalent: values must be positive, got max = $(v̄)")
    acc = zero(float(v̄))
    @inbounds for i in eachindex(values, weights)
        acc += weights[i] * (values[i]/v̄)^(1-σ)
    end
    return v̄ * acc^(1/(1-σ))
end

function create_initial_guess()
    function init_vDefault(state)
        y = state[end]
        return h(y)
    end

    function init_vRepay(state)
        b, y = state[1], state[end]
        return y-b
    end

    function init_q()
        return 1.0
    end

    function init_pol()
        return 0.0
    end

    polASG = AdaptiveSparseGrid(inputs -> init_pol(), 
        [b_L, y_L], [b_H, y_H], max_depth = 7, tol = 1e-8)

    qASG = AdaptiveSparseGrid(inputs -> init_q(), 
        [b_L, y_L], [b_H, y_H], max_depth = 7, tol = 1e-8)

    vDefaultASG = AdaptiveSparseGrid(y -> init_vDefault(y), 
        [y_L], [y_H], max_depth = 8, tol = 1e-8)

    vRepayASG = AdaptiveSparseGrid(state -> init_vRepay(state), 
        [b_L, y_L], [b_H, y_H], max_depth = 7, tol = 1e-8)

    return polASG, qASG, vDefaultASG, vRepayASG
end
