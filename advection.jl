using LinearAlgebra
using DelimitedFiles
using Statistics
using CairoMakie

# physical properties
ρ = 1000 # water density [kg/m^3]
μ = 1.0e-3 # dynamic viscosity [Pa s]
g = 9.81 # gravitational acceleration [m/s^2]

# model domain
Lx = 100 # domain length in x direction [m]
Ly = 50 # domain length in y direction [m]

nx = 100 # number of cells in x direction
ny = 50 # number of cells in y direction
nc = nx * ny # total number of cells

dx = Lx / nx # = 1 m
dy = Ly / ny # = 1 m

# boundary conditions
uin = 0.1 / 100 / 60 # [m/s] == 0.1 cm/min <-- Neumann on left boundary
hR = 9 # [m] <-- Dirichlet on right boundary

## transmissibility matrix
# read in permeability field from CSV file
kmat = readdlm("data/permeability_field.csv", ',', Float64) # [m^2]
kmean = mean(log10.(reshape(kmat, nc)))
kstd = std(log10.(reshape(kmat, nc)))

Y = (log10.(kmat) .- kmean) ./ kstd # standardized spatial field
σ = 0.1 # heterogeneity level
log10k = kmean .+ σ .* Y # rescaled log10(k) field
Kmat = ρ*g/μ .* 10 .^ log10k # hydraulic conductivity [m/s]

K_har_x = (1 ./ Kmat[:, 1:(nx - 1)] .+ 1 ./ Kmat[:, 2:nx]) .^ -
1 # harmonic mean in x direction
K_har_y = (1 ./ Kmat[1:(ny - 1), :] .+ 1 ./ Kmat[2:ny, :]) .^ -
1 # harmonic mean in y direction

# TPFA
TLmat = fill(NaN, ny, nx)
TLmat[:, 1] .= 0.0 # left boundary (Neumann): flux enters through b
TLmat[:, 2:nx] .= (2*dy/dx) .* K_har_x # transmissibility in x direction

TRmat = fill(NaN, ny, nx)
TRmat[:, nx] .= (2*dy/dx) .* Kmat[:, nx] # right boundary (Dirichlet)
TRmat[:, 1:(nx - 1)] .= (2*dy/dx) .* K_har_x # transmissibility in x direction

TUmat = fill(NaN, ny, nx)
TUmat[1, :] .= 0.0 # top boundary: no-flow
TUmat[2:ny, :] .= (2*dx/dy) .* K_har_y # transmissibility in y direction

TDmat = fill(NaN, ny, nx)
TDmat[ny, :] .= 0.0 # bottom boundary: no-flow
TDmat[1:(ny - 1), :] .= (2*dx/dy) .* K_har_y # transmissibility in y direction

TL = reshape(TLmat, nc)
TR = reshape(TRmat, nc)
TU = reshape(TUmat, nc)
TD = reshape(TDmat, nc)

TC = TL .+ TR .+ TU .+ TD # transmissibility matrix diagonal

Tmat =
    diagm(0 => TC) +
    diagm(-ny => -TL[(ny + 1):nc]) +
    diagm(-1 => -TU[2:nc]) +
    diagm(1 => -TD[1:(nc - 1)]) +
    diagm(ny => -TR[1:(nc - ny)])

## RHS vector b
bmat = zeros(ny, nx)
bmat[:, 1] .= uin * dy # left boundary: specified flux (Neumann)
bmat[:, nx] .= TRmat[:, nx] .* hR # right boundary: specified head (Dirichlet)
b = reshape(bmat, nc)

## head vector h
h = Tmat \ b
hmat = reshape(h, ny, nx)

## Darcy velocity
dh_x = hmat[:, 1:(nx - 1)] .- hmat[:, 2:nx] # head difference in x direction
dh_y = hmat[1:(ny - 1), :] .- hmat[2:ny, :] # head difference in y direction

UL = fill(NaN, ny, nx)
UL[:, 1] .= uin # left boundary: specified inflow
UL[:, 2:nx] .= TLmat[:, 2:nx] ./ dy .* dh_x

UR = fill(NaN, ny, nx)
UR[:, nx] .= TRmat[:, nx] ./ dy .* (hmat[:, nx] .- hR) # right boundary: outflow to hR
UR[:, 1:(nx - 1)] .= TRmat[:, 1:(nx - 1)] ./ dy .* dh_x

UU = fill(NaN, ny, nx)
UU[1, :] .= 0.0 # top boundary: no-flow
UU[2:ny, :] .= TUmat[2:ny, :] ./ dx .* dh_y

UD = fill(NaN, ny, nx)
UD[ny, :] .= 0.0 # bottom boundary: no-flow
UD[1:(ny - 1), :] .= TDmat[1:(ny - 1), :] ./ dx .* dh_y

Ux = 0.5 * (UL .+ UR) # cell-centered Darcy velocity, x direction
Uy = 0.5 * (UU .+ UD) # cell-centered Darcy velocity, y direction

## visualization
set_theme!(theme_latexfonts())
update_theme!(
    fontsize = 18,
    Axis = (yreversed = true, aspect = DataAspect()), # yreversed mimics imagesc
)
X = range(dx / 2, step = dx, length = nx) # cell-center coordinates
Ycoord = range(dy / 2, step = dy, length = ny) # `Y` is the standardized field,
don't overwrite it
skip = 5

# Darcy velocities are tiny relative to the domain, so the raw vectors are
# invisible; lengthscale blows them up to a few cells long
maxU = maximum(sqrt.(Ux .^ 2 .+ Uy .^ 2))
add_velocity!(ax) = arrows2d!(
    ax,
    X[1:skip:end],
    Ycoord[1:skip:end],
    Ux[1:skip:end, 1:skip:end]',
    Uy[1:skip:end, 1:skip:end]',
    color = :black,
    lengthscale = 2 * skip / maxU,
)

fig = Figure(size = (800, 800))
ax1 = Axis(fig[1, 1], title = "log10(K)")
hm1 = heatmap!(ax1, X, Ycoord, log10.(Kmat)')
Colorbar(fig[1, 2], hm1)

ax2 = Axis(fig[2, 1], xlabel = "x [m]", ylabel = "y [m]", title = "Hydraulic head")
hm2 = heatmap!(ax2, X, Ycoord, hmat')
Colorbar(fig[2, 2], hm2, label = "[m]")
add_velocity!(ax2)
display(fig)

## advection
ϕ = 0.3 # porosity
day = 24 * 60 * 60 # [s]
dt_max = 0.1 * min(dx, dy) / maximum(sqrt.(Ux .^ 2 .+
Uy .^ 2) ./ ϕ) # CFL on pore velocity

# one explicit upwind step; returns the new concentration field
function advect(CC, C_in, dt)
    CL = [fill(C_in, ny, 1) CC[:, 1:(nx - 1)]] # left boundary condition of C_in
    CR = [CC[:, 2:nx] CC[:, nx]] # zero-
    gradient concentration condition on outflow boundary
    CU = [CC[1:1, :]; CC[1:(ny - 1), :]]
    CD = [CC[2:ny, :]; CC[ny:ny, :]]

    # upwind: take the concentration from the side the flow comes from
    Fadv_L = ifelse.(UL .>= 0, UL .* CL, UL .* CC)
    Fadv_U = ifelse.(UU .>= 0, UU .* CU, UU .* CC)
    Fadv_D = ifelse.(UD .>= 0, UD .* CC, UD .* CD)
    Fadv_R = ifelse.(UR .>= 0, UR .* CC, UR .* CR)

    netFadv = -Fadv_L .* dy .- Fadv_U .* dx .+ Fadv_D .* dx .+ Fadv_R .* dy
    return CC .- dt / ϕ / dx / dy .* netFadv
end

# animation figure: heatmaps read from Observables so record can update them
CCobs = Observable(zeros(nx, ny))
titleobs = Observable("time: 0 days")
figC = Figure(size = (800, 800))
axK = Axis(figC[1, 1], title = "log10(K)")
hmK = heatmap!(axK, X, Ycoord, log10.(Kmat)')
Colorbar(figC[1, 2], hmK)
add_velocity!(axK)
axC = Axis(figC[2, 1], xlabel = "x [m]", ylabel = "y [m]", title = titleobs)
hmC = heatmap!(axC, X, Ycoord, CCobs, colorrange = (0, 1))
Colorbar(figC[2, 2], hmC)

# `let` keeps t and CC local so the do-block can update them
CC = let t = 0.0, CC = zeros(ny, nx)
    record(figC, "figures/advection.mp4", 1:30; framerate = 5) do iday
        # step up to the end of this day, shortening the last step to land on it
        while t < iday * day
            C_in = t < 10 * day ? 1.0 : 0.0 # stop injecting after 10 days
            dt = min(dt_max, iday * day - t)
            CC = advect(CC, C_in, dt)
            t += dt
        end
        CCobs[] = CC'
        titleobs[] = "time: $iday days"
    end
    CC
end
