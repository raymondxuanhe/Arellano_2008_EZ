pref_label() = "EIS = $(round(1/ρ, digits=2)), RA = $(round(σ, digits=2))"

function plot_prices(y, qASG)
    q(b) = qASG([b, y])
    plot(bgrid, q.(bgrid), xlabel = "B′", ylabel = "Price",
        title = "Pricing Functions for y = $(round(y, digits=2)) ($(pref_label()))",
        dpi = 300, grid = false, label = false)
end

function plot_policy(polASG)
    pol_ss(b) = polASG([b, 1.0])
    pol_l(b) = polASG([b, 0.95])
    pol_h(b) = polASG([b, 1.05])

    plot(bgrid, pol_ss.(bgrid), label="y=1.0", xlabel="b", ylabel="b′",
        title = "Policy function ($(pref_label()))", dpi = 300)
    plot!(bgrid,pol_l.(bgrid),label="y=0.95")
    plot!(bgrid,pol_h.(bgrid),label="y=1.05")
end

function plot_vDefault_util(vDefaultASG)
    vD(y) = util(vDefaultASG([y]))
    plot(ygrid, vD.(ygrid), xlabel = "y", ylabel = "Value (Util Units)",
        title = "vDefault ($(pref_label()))", label = false, dpi = 300, grid = false)
end

function plot_vDefault(vDefaultASG)
    vD(y) = vDefaultASG([y])
    plot(ygrid, vD.(ygrid), xlabel = "y", ylabel = "Value (CE Units)",
        title = "vDefault ($(pref_label()))", label = false, dpi = 300, grid = false)
end

function plot_vRepay_util(vRepayASG)
    vR_low(y) = util(vRepayASG([0.0, y]))
    vR_med(y) = util(vRepayASG([0.2, y]))
    vR_high(y) = util(vRepayASG([0.4, y]))

    plot(ygrid, vR_low.(ygrid), label = "vRepay (No Debt)", xlabel = "y",
        ylabel = "Value (Util Units)", title = "vRepayment ($(pref_label()))",
        legend = :bottomright, dpi = 300, grid = false)
    plot!(ygrid, vR_med.(ygrid), label = "vRepay (Medium Debt)")
    plot!(ygrid, vR_high.(ygrid), label = "vRepay (High Debt)")
end

function plot_vRepay(vRepayASG)
    vR_low(y) = vRepayASG([0.0, y])
    vR_med(y) = vRepayASG([0.2, y])
    vR_high(y) = vRepayASG([0.4, y])

    plot(ygrid, vR_low.(ygrid), label = "vRepay (No Debt)", xlabel = "y",
        ylabel = "Value (CE Units)", title = "vRepayment ($(pref_label()))",
        legend = :bottomright, dpi = 300, grid = false)
    plot!(ygrid, vR_med.(ygrid), label = "vRepay (Medium Debt)")
    plot!(ygrid, vR_high.(ygrid), label = "vRepay (High Debt)")
end
