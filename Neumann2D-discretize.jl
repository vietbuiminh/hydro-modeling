using LinearAlgebra
using Plots

# physical properties
ρ = 1000 # water density [g/cm^3]
μ = 1.0e-3 # dynamic viscosity [Pa s]
g = 9.81 # gravitational acceleration [m/s^2]
k = 1e-10 # instrinsic permeability [m^2]

# model domain
Lx = 100 # domain length in x direction [m]
Ly = 50 # domain length in y direction [m]

nx = 100 # number of cells in x direction
ny = 50 # number of cells in y direction
nc = nx * ny # total number of cells

dx = Lx / nx
dy = Ly / ny

# boundary conditions
u_in = 0.1 / 100 / 60 # left boundary flux [m/s] = 0.1 cm/min [Neumann type]
hR = 9 # right boundary head [m] [Dirichlet type]

# transmissibility matrix
# in class the code contain an input file for Kmat, but here we just use a constant
# value for Kmat
# TODO: read in Kmat from a file, and then calculate the harmonic mean in x and y
# directions
kmat = readmatrix("data/Kmat.csv") # read in Kmat from a file
# Kmat = ρ*g*k/μ .* ones(ny, nx) # m/s <--- this is K * matrix
Kmat = ρ*g*kmat/μ .* kmat

kmean = mean(log10(reshape(kmat, nc)))
kstd = std(log10(reshape(kmat, nc)))

Y = (log10(reshape(kmat, nc)) - kmean) / kstd
kmat = kmean .+ 0.1 .* Y # 0.1 or 1.0 and make some comment
Kmat = ρ*g*kmat/μ .* 10 .^ kmat

K_har_x = (1 ./ Kmat[:, 1:(nx - 1)] .+ 1 ./ Kmat[:, 2:nx]) .^ -
1 # harmonic mean in x direction
K_har_y = (1 ./ Kmat[1:(ny - 1), :] .+ 1 ./ Kmat[2:ny, :]) .^ -
1 # harmonic mean in y direction

TLmat = fill(NaN, ny, nx)
TLmat[:, 1] .= 0.0 # left boundary: no transmissibility
TLmat[:, 2:nx] .= (2/dx*dy) .* K_har_x # transmissibility in x direction

TRmat = fill(NaN, ny, nx)
TRmat[:, nx] .= (2/dx*dy) .* Kmat[:, nx]
TRmat[:, 1:(nx - 1)] .= (2/dx*dy) .* K_har_x # transmissibility in x direction

TUmat = fill(NaN, ny, nx)
TUmat[1, :] .= 0.0 # top boundary: no transmissibility
TUmat[2:ny, :] .= (2/dx*dy) .* K_har_y # transmissibility in y direction

TDmat = fill(NaN, ny, nx)
TDmat[ny, :] .= 0.0 # bottom boundary: no transmissibility
TDmat[1:(ny - 1), :] .= (2/dx*dy) .* K_har_y # transmissibility in y direction

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
bmat[:, 1] .= u_in .* dy # left boundary: specified flux (Neumann)
bmat[:, nx] .= TRmat[:, nx] .* hR # right boundary: specified head (Dirichlet)

b = reshape(bmat, nc) # RHS vector

# hydraulic head vector h
h = Tmat \ b # solve for hydraulic head

hmat = reshape(h, ny, nx) # reshape h vector to matrix form
heatmap(hmat, xlabel = "x cell", ylabel = "y cell", title = "Hydraulic head")

# TODO make the subplot of 2 where the first show the kmat and the bottom is the result
# of hmat
p = plot(layout = (2, 1))
heatmap!(
    p[1],
    log10(Kmat),
    xlabel = "x cell",
    ylabel = "y cell",
    title = "Log10(Permeability)",
)
heatmap!(p[2], hmat, xlabel = "x cell", ylabel = "y cell", title = "Hydraulic head")

# Darcy velocity 

dh_x = hmat[:, 1:(nx - 1)] .- hmat[:, 2:nx] # head difference in x direction
dh_y = hmat[1:(ny - 1), :] .- hmat[2:ny, :] # head difference in y direction

UL = fill(NaN, ny, nx)
UL[:, 1] .= u_in # left boundary: specified flux (Neumann)
UL[:, 2:nx] .= TLmat[:, 2:nx] .* dh_x ./ dx # Darcy velocity in x direction

UR = fill(NaN, ny, nx)
UR[:, nx] .= TRmat[:, nx] .* dh_x[:, nx -
1] ./ dx # right boundary: specified head (Dirichlet)
UR[:, 1:(nx - 1)] .= TRmat[:, 1:(nx - 1)] .* dh_x ./ dx # Darcy velocity in x direction

UU = fill(NaN, ny, nx)
UU[1, :] .= 0.0 # top boundary: no transmissibility
UU[2:ny, :] .= TUmat[2:ny, :] .* dh_y ./ dx # Darcy velocity in y direction

UD = fill(NaN, ny, nx)
UD[ny, :] .= 0.0 # bottom boundary: no transmissibility
UD[1:(ny - 1), :] .= TDmat[1:(ny - 1), :] .* dh_y ./ dx # Darcy velocity in y direction

Ux = 0.5 * (UL .+ UR) # average Darcy velocity in x direction
Uy = 0.5 * (UU .+ UD) # average Darcy velocity in y direction
