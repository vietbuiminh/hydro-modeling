using Statistics
using Printf
using CairoMakie

include("funcs/flow_transport_2d.jl")

# HW5 Problem 1: advection-dispersion in a homogeneous aquifer

# physical properties
ρ = 1000 # water density [kg/m^3]
μ = 1.0e-3 # dynamic viscosity [Pa s]
g = 9.81 # gravitational acceleration [m/s^2]
k = 1.0e-10 # intrinsic permeability [m^2]

ϕ = 0.3 # porosity [-]
αVG = 1.0 # inverse air entry scale [1/m]
nVG = 2 # pore size distribution index [-]
mVG = 1 - 1 / nVG # van Genuchten parameter [-]
lVG = 0.5 # pore connectivity parameter [-]
θs = ϕ # saturated water content [-]
θr = 0.05 # irreducible residual water content [-]

D_diff = 1.0e-9 # molecular diffusion coeff. [m^2/s]
α = 1.0 # isotropic dispersivity [m]

# model domain
Lx = 100 # domain length in x direction [m]
Ly = 50 # domain length in y direction [m]
nx = 100 # number of cells in x direction
ny = 50 # number of cells in y direction
dx = Lx / nx # = 1 m
dy = Ly / ny # = 1 m

# boundary conditions
# hL = 10 # left boundary head [m] (Dirichlet)
# hR = 9 # right boundary head [m] (Dirichlet)
# day = 24 * 60 * 60 # [s]
# c_in(t) = t < 10 * day ? 1.0 : 0.0 # inject c = 1 for the first 10 days
uin = 0.1/100/60 # [m/s] = 0.1 cm/min <-- Neumann boundary top boundary
psiB = 0 # bottom boundary head [m] (Dirichlet) implies water table at the bottom of the domain

## flow
# Kmat = (ρ*g*k/μ) .* ones(ny, nx) # hydraulic conductivity [m/s]
kmat = readdlm("data/permeability_field.csv", ',', Float64) # [m^2]
kmean = mean(log10.(kmat))
kstd = std(log10.(kmat))

Y = (log10.(kmat) .- kmean) ./ kstd
kmat = kmean + 0.1*Y

KmatSat = (ρ*g*kmat/μ) .* 10 .^ kmat # saturated hydraulic conductivity [m/s]

# initial settings
t = 0 # initial time [s]
tEnd = 3*24*60*60 # 3 days [s]
