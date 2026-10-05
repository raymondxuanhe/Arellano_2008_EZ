# Arellano (2008) with Epstein-Zin Preferences

Julia code solving the Arellano (2008) sovereign default model, extended to allow for Epstein-Zin preferences instead of standard expected utility (CRRA).

Separating risk aversion from the elasticity of intertemporal substitution (EIS) typically makes these models numerically delicate: the recursion over the continuation value and bond price can become unstable as risk aversion and the EIS diverge from each other. This code is built to remain numerically robust across a wide range of EIS and risk aversion parameter combinations.

## Structure

- `arellano_2008_EZ.jl` — main driver script: sets preference parameters, builds the state space, runs the value function iteration (`solve_model`), and plots results.
- `arellano_2008_EZ_parameters.jl` — model and grid parameters (income process, bond grid, discretization).
- `arellano_2008_EZ_auxiliary_functions.jl` — helper functions used across the solution routines.
- `arellano_2008_EZ_updating_functions.jl` — Bellman updating steps for bond prices, default/repayment value functions, and policy functions.
- `arellano_2008_EZ_plotting_functions.jl` — plotting utilities for prices, policies, and value functions.
- `arellano_2008_EZ_validation_functions.jl` — validation/diagnostic routines.

## Method

The model is solved via value function iteration using adaptive sparse grids ([AdaptiveSparseGrids.jl](https://github.com/)) over the (debt, income) state space, with Gauss-Hermite quadrature for integrating over the income shock.

## Requirements

Dependencies are listed in `Project.toml`. Instantiate the environment with:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Usage

```julia
include("arellano_2008_EZ.jl")
```

This solves the model and produces plots of bond prices, the default/repayment policy, and the value functions.
