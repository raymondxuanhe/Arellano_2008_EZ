#Solving the Arellano (2008) model of sovereign default with Epstein-Zin preferences
#Author: Raymond He
#Last updated: October 5, 2026

using AdaptiveSparseGrids, FastGaussQuadrature, LinearAlgebra, StaticArrays
using Optimization, OptimizationNLopt, OptimizationMultistartOptimization, SciMLBase
using Parameters, Plots, Printf
plotlyjs()

const D = 2;

@with_kw mutable struct Validation
    polVal::Array{Float64, D}
    qVal::Array{Float64, D}
    vDVal::Array{Float64, 1}
    vRVal::Array{Float64, D}
end

include("arellano_2008_EZ_parameters.jl");
include("arellano_2008_EZ_auxiliary_functions.jl");
include("arellano_2008_EZ_updating_functions.jl");
include("arellano_2008_EZ_plotting_functions.jl");
include("arellano_2008_EZ_validation_functions.jl");

const σ = 30.0;
const ρ = 1.1;

const α = 0.02;

const states = [SVector(b, y) for b in bgrid, y in ygrid];

function solve_model(max_iter, tol)
    #Create initial guess
    iiter = 0
    polASG, qASG, vDefaultASG, vRepayASG = create_initial_guess();

    ValA = create_init_val(states, polASG, qASG, vDefaultASG, vRepayASG);

    while iiter < max_iter
        iiter += 1

        qASG_new = AdaptiveSparseGrid(inputs -> update_q(inputs, vDefaultASG, vRepayASG), 
            [b_L, y_L], [b_H, y_H], max_depth = 7, tol = 1e-8)

        vDefaultASG_new = AdaptiveSparseGrid(inputs -> update_vDefault(inputs, vDefaultASG, vRepayASG), 
            [y_L], [y_H], max_depth = 7, tol = 1e-8)

        RASG = AdaptiveSparseGrid(inputs -> expected_CV(inputs, vDefaultASG, vRepayASG),
            [b_L, y_L], [b_H, y_H], max_depth = 7, tol = 1e-8)

        vRepayASG_new = AdaptiveSparseGrid(inputs -> update_vRepay(inputs, qASG_new, RASG),
            [b_L, y_L], [b_H, y_H], max_depth = 7, tol = 1e-8)

        polASG_new = AdaptiveSparseGrid(inputs -> update_pol(inputs, qASG_new, RASG),
            [b_L, y_L], [b_H, y_H], max_depth = 7, tol = 1e-8)

        #Check for convergence
        dp, dq, dvD, dvR = update_and_diff!(ValA, states, polASG_new, qASG_new, vDefaultASG_new, vRepayASG_new)

        @printf("%5d %20.5e %20.5e %20.5e %20.5e\n", iiter, dp, dq, dvD, dvR)

        polASG = deepcopy(polASG_new)
        qASG = deepcopy(qASG_new)
        vDefaultASG = deepcopy(vDefaultASG_new)
        vRepayASG = deepcopy(vRepayASG_new)

        if dvR < tol
            println("The model converged after $iiter iterations!")
            break
        end
    end
    return polASG, qASG, vDefaultASG, vRepayASG
end

polASG, qASG, vDefaultASG, vRepayASG = solve_model(100, 1e-5)

plot_prices(0.9, qASG)
plot_prices(1.0, qASG)
plot_prices(1.1, qASG)

plot_policy(polASG)

plot_vDefault_util(vDefaultASG)
plot_vRepay_util(vRepayASG)
