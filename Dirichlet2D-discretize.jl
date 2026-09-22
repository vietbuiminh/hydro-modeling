using LinearAlgebra
using CairoMakie

# HW4 Problem 1-2: homogeneous k, Dirichlet BCs at left and right,
# no-flow (Neumann, zero-flux) BCs at top and bottom

# physical properties
ρ = 1000 # water density [kg/m^3]
μ = 1.0e-3 # dynamic viscosity [Pa s]
g = 9.81 # gravitational acceleration [m/s^2]
k = 1e-10 # intrinsic permeability [m^2]

# model domain
Lx = 100 # domain length in x direction [m]
Ly = 50 # domain length in y direction [m]

nx = 6 # number of cells in x direction
ny = 3 # number of cells in y direction
nc = nx * ny # or N in note total number of cells = 18

dx = Lx / nx # = 50/3 m
dy = Ly / ny # = 50/3 m

# boundary conditions
hL = 10 # left boundary head [m] [Dirichlet type]
hR = 9 # right boundary head [m] [Dirichlet type]

# hydraulic conductivity (homogeneous)
Kmat = (ρ*g*k/μ) .* ones(ny, nx) # m/s

K_har_x = (1 ./ Kmat[:, 1:(nx - 1)] .+ 1 ./ Kmat[:, 2:nx]) .^ -
1 # harmonic mean in x direction
K_har_y = (1 ./ Kmat[1:(ny - 1), :] .+ 1 ./ Kmat[2:ny, :]) .^ -
1 # harmonic mean in y direction

# transmissibility matrix
# TPFA
TLmat = fill(NaN, ny, nx)
TLmat[:, 1] .= (2*dy/dx) .* Kmat[:, 1] # left boundary (Dirichlet)
TLmat[:, 2:nx] .= (2*dy/dx) .* K_har_x # transmissibility in x direction

TRmat = fill(NaN, ny, nx)
TRmat[:, nx] .= (2*dy/dx) .* Kmat[:, nx] # right boundary (Dirichlet)
TRmat[:, 1:(nx - 1)] .= (2*dy/dx) .* K_har_x # transmissibility in x direction

TUmat = fill(NaN, ny, nx)
TUmat[1, :] .= 0.0 # top boundary: no transmissibility
TUmat[2:ny, :] .= (2*dy/dx) .* K_har_y # transmissibility in y direction

TDmat = fill(NaN, ny, nx)
TDmat[ny, :] .= 0.0 # bottom boundary: no transmissibility
TDmat[1:(ny - 1), :] .= (2*dy/dx) .* K_har_y # transmissibility in y direction

TL = reshape(TLmat, nc)
TR = reshape(TRmat, nc)
TU = reshape(TUmat, nc)
TD = reshape(TDmat, nc)

TC = TL .+ TR .+ TU .+ TD # transmissibility matrix diagonal

Tmat =
    diagm(0 => TC) +
    diagm(-1 => -TU[2:nc]) +
    diagm(1 => -TD[1:(nc - 1)]) +
    diagm(-ny => -TL[(ny + 1):nc]) +
    diagm(ny => -TR[1:(nc - ny)]) # m/s

# RHS vector b
bmat = zeros(ny, nx)
bmat[:, 1] .= TLmat[:, 1] .* hL # left boundary: specified head (Dirichlet)
bmat[:, nx] .= TRmat[:, nx] .* hR # right boundary: specified head (Dirichlet)
b = reshape(bmat, nc) # RHS vector

# hydraulic head vector h (Problem 2)
h = Tmat \ b # solve for hydraulic head

hmat = reshape(h, ny, nx)

# plot
set_theme!(theme_latexfonts())
update_theme!(
    fontsize = 20,
    figure_padding = (10, 50, 10, 10),
    Axis = (yreversed = true,),
)

fig = Figure(size = (700, 350))
ax = Axis(
    fig[1, 1],
    xlabel = "Lx (m)",
    ylabel = "Ly (m)",
    title = "Hydraulic head",
)
X = range(dx / 2, step = dx, length = nx)
Y = range(dy / 2, step = dy, length = ny)
hm = heatmap!(ax, X, Y, hmat')
Colorbar(fig[1, 2], hm, label = "(m)")

display(fig)
save("figures/hydraulic_head.png", fig)
