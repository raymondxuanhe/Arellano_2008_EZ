function update_q(inputs, vDefaultASG, vRepayASG)
    b′, y = inputs[1], inputs[end]
    q = 0.0
    for i in 1:nquad
        y′ = clamp(exp(ρy*log(y) + y_nodes[i]), y_L, y_H)
        vD = vDefaultASG([y′])
        vR = vRepayASG([b′, y′])
        uD, uR = util(vD), util(vR)
        m = max(uD, uR)
        repay_prob = exp((uR-m)/α)/(exp((uR-m)/α) + exp((uD-m)/α))
        q += y_weights[i] * repay_prob
    end
    return clamp(q/(1+r), q_L, q_H)
end

function expected_CV(inputs, vDefaultASG, vRepayASG)
    b′, y = inputs[1], inputs[end]

    values  = Vector{Float64}(undef, nquad)
    weights = Vector{Float64}(undef, nquad)
    for j in 1:nquad
        y′ = clamp(exp(ρy*log(y) + y_nodes[j]), y_L, y_H)
        vD = vDefaultASG([y′])
        vR = vRepayASG([b′, y′])
        W  = util_inv(ex_ante_value(util(vD), util(vR)))
        values[j]  = W
        weights[j] = y_weights[j]
    end
    return certainty_equivalent(values, weights)
end

function update_vDefault(point, vDefaultASG, vRepayASG)
    y = point[end]
    flow_term = h(y)^(1-ρ)

    values  = Vector{Float64}(undef, 2*nquad)
    weights = Vector{Float64}(undef, 2*nquad)
    for j in 1:nquad
        y′ = clamp(exp(ρy*log(y) + y_nodes[j]), y_L, y_H)
        vD = vDefaultASG([y′])
        vR = vRepayASG([0.0, y′])
        W  = util_inv(ex_ante_value(util(vD), util(vR)))

        values[j]           = vD
        weights[j]           = y_weights[j] * (1-θ)
        values[nquad + j]   = W
        weights[nquad + j]   = y_weights[j] * θ
    end

    R = certainty_equivalent(values, weights)
    return ((1-β)*flow_term + β*R^(1-ρ))^(1/(1-ρ))
end

function update_vRepay(point, qASG_new, RASG)
    function obj_fn(x, state)
        b, y = state[1], state[end]
        b′ = x[1]

        q = clamp(qASG_new([b′, y]), q_L, q_H)

        c = y - b + q*b′
        flow_term = c^(1-ρ)

        R = RASG([b′, y])
        CV_term = R^(1-ρ)

        return -(((1-β)*flow_term + β*CV_term)^(1/(1-ρ)))
    end
    x0 = [point[1]/2.0]
    prob = SciMLBase.OptimizationProblem(obj_fn, x0,
            point, lb = [b_L], ub = [b_H])
    sol = solve(prob, MultistartOptimization.TikTak(25), NLopt.LN_NELDERMEAD())
    return -sol.objective
end

function update_pol(point, qASG_new, RASG)
    function obj_fn(x, state)
        b, y = state[1], state[end]
        b′ = x[1]

        q = clamp(qASG_new([b′, y]), q_L, q_H)

        c = y - b + q*b′
        flow_term = c^(1-ρ)

        R = RASG([b′, y])
        CV_term = R^(1-ρ)

        return -(((1-β)*flow_term + β*CV_term)^(1/(1-ρ)))
    end
    x0 = [point[1]/2.0]
    prob = SciMLBase.OptimizationProblem(obj_fn, x0,
            point, lb = [b_L], ub = [b_H])
    sol = solve(prob, MultistartOptimization.TikTak(25), NLopt.LN_NELDERMEAD())
    return sol.u[1]
end
