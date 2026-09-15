using LinearAlgebra
using Plots

# physical properties
ρ = 1000 # g/cm^3
g = 9.81 # m/s^2
μ = 1.0e-3 # Pa*s

k = 1e-12 # m^2

# model dimension
Lx = 100 # domain dimention length in x direction [m]
Ly = 50 # domain dimention length in y direction [m]

nx = 6 # number of cells in x direction
ny = 3 # number of cells in y direction

dx = Lx / nx
dy = Ly / ny

# boundaries conditions
uL = 0.1 / 100 / 60 # left boundary flux [m/s] = 0.1 cm/min [Neumann type]
hR = 1 # head value on the right boundary [m] [Dirichlet type]

# transmissibility matrix
K = ρ*g*k/μ # m/s
Kmat = K*ones(ny, nx) # m/s

TL = fill(NaN, ny, nx)
TL[:, 1] .= 0.0 # m/s
TL[:, 2:nx] .= 2*dy/dx * (1 ./ Kmat[:, 1:(nx - 1)] + 1 ./ Kmat[:, 2:nx]) .^ -1 # m/s

TR = fill(NaN, ny, nx)
TR[:, nx] .= 2*dy/dx*Kmat[:, nx] # m/s
TR[:, 1:(nx - 1)] .= 2*dy/dx * (1 ./ Kmat[:, 1:(nx - 1)] + 1 ./ Kmat[:, 2:nx]) .^ -
1 # m/s

TU = fill(NaN, ny, nx)
TU[1, :] .= 0.0
TU[2:ny, :] .= 2*dx/dy * (1 ./ Kmat[1:(ny - 1), :] + 1 ./ Kmat[2:ny, :]) .^ -1 # m/s   

TD = fill(NaN, ny, nx)
TD[ny, :] .= 0
TD[1:(ny - 1), :] .= 2*dx/dy * (1 ./ Kmat[1:(ny - 1), :] + 1 ./ Kmat[2:ny, :]) .^ -
1 # m/s

TC = TL .+ TR .+ TU .+ TD # m/s
TC = reshape(TC, nx * ny)

TL = reshape(TL, nx*ny)
TU = reshape(TU, nx*ny)
TR = reshape(TR, nx*ny)
TD = reshape(TD, nx*ny)

Tmat =
    diagm(0 => TC) +
    diagm(-1 => -TU[2:(nx * ny)]) +
    diagm(1 => -TD[1:(nx * ny - 1)]) +
    diagm(-ny => -TL[(ny + 1):(nx * ny)]) +
    diagm(ny => -TR[1:(nx * ny - ny)]) # m/s

# RHS vector:b
b = zeros(ny, nx)
b[:, 1] .= uL*dy
b[:, nx] .= TR[((nx - 1) * ny + 1):(nx * ny)] * hR
b = reshape(b, nx*ny)

TL = reshape(TL, nx*ny)
TR = reshape(TR, nx*ny)
TU = reshape(TU, nx*ny)
TD = reshape(TD, nx*ny)
TC = reshape(TC, nx*ny)

# PLOT HEATMAP
plot(
    reshape(TC, nx, ny),
    seriestype = :heatmap,
    xlabel = "x",
    ylabel = "y",
    title = "Transmissibility",
)
