using LinearAlgebra
using Plots

# physical properties
ρ = 1000 # g/cm^3
g = 9.81 # m/s^2
μ = 1.0e-3 # Pa*s

k = 1e-12 # m^2

# model dimension
Lx = 1
nx = 4
dx = Lx / nx

# boundaries conditions
hL = 10
hR = 9

# transmissibility
K = ρ*g*k/μ # m/s
Kmat = fill(K, nx, 1) # m/s

TL = fill(NaN, nx)
TL[1] = 2/dx*Kmat[1] # m/s
TL[2:nx] = 2/dx * (1 ./ Kmat[1:(nx - 1)] + 1 ./ Kmat[2:nx]) .^ -1 # m/s

TR = fill(NaN, nx)
TR[nx] = 2/dx*Kmat[nx] # m/s
TR[1:(nx - 1)] = 2/dx * (1 ./ Kmat[1:(nx - 1)] + 1 ./ Kmat[2:nx]) .^ -1 # m/s

TC = TL + TR # m/s

Tmat = Matrix(Tridiagonal(-TL[2:nx], TC, -TR[1:(nx - 1)])) # m/s

# RHS vector:b 
b = zeros(nx)
b[1] = TL[1]*hL
b[nx] = TR[nx]*hR

# estimate h = T/b
h = Tmat \ b
xcell = range(dx/2, stop = Lx - dx/2, length = nx)

Plots.plot(
    collect(xcell),
    h,
    seriestype = :scatter,
    label = "h",
    xlabel = "x (m)",
    ylabel = "h (m)",
    title = "1D Discretization of Flow in Porous Media",
    legend = :topright,
)
