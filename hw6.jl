using Statistics
using Printf
using CairoMakie

include("funcs/variable_density_2d.jl")

# HW6: variable-density flow and solute transport in a homogeneous aquifer

# physical properties
ρ = 1000 # freshwater density [kg/m^3]
μ = 1.0e-3 # dynamic viscosity [Pa s]
g = 9.81 # gravitational acceleration [m/s^2]
k = 1.0e-10 # intrinsic permeability [m^2]

ϕ = 0.3 # porosity [-]
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
uin = 0.1 / 100 / 60 # left boundary inflow [m/s] == 0.1 cm/min (Neumann)
hR = 9 # right boundary freshwater head [m] (Dirichlet)
day = 24 * 60 * 60 # [s]
c_in(t) = t < 10 * day ? 1.0 : 0.0 # inject c = 1 for the first 10 days

Kmat = (ρ*g*k/μ) .* ones(ny, nx) # hydraulic conductivity [m/s]
K = Kmat[1]

## three density contrasts, all evaluated at day 7
Δρ_ls = [0, 15, 30] # [kg/m^3]
β_ls = Δρ_ls ./ ρ # relative density contrast [-]
t7 = 7 * day
runs = [solve_vd(Kmat, dx, dy, uin, hR, β, ϕ, D_diff, α, c_in, [t7])[1] for β in β_ls]

X = range(dx / 2, step = dx, length = nx) # cell-center coordinates [m]
Ycoord = range(dy / 2, step = dy, length = ny) # y = 0 at the top, +y downward

## Problem 2: maximum velocities at day 7 (cell-centered Ux, Uy)
maxUy = [maximum(abs, r.Uy) for r in runs]
maxU = [maximum(sqrt.(r.Ux .^ 2 .+ r.Uy .^ 2)) for r in runs]

println("| Δρ (kg/m³) | β (-) | max|Uy| (m/s) | max|U| (m/s) |")
println("|---|---|---|---|")
for (Δρ, β, uy, u) in zip(Δρ_ls, β_ls, maxUy, maxU)
    @printf("| %d | %.3f | %.3e | %.3e |\n", Δρ, β, uy, u)
end
@printf("max|Uy|(30) / max|Uy|(15) = %.3f (linear would be 2)\n", maxUy[3] / maxUy[2])

## Problem 1: horizontal transport distance in the top and bottom rows
# x position of the c = 0.5 front in a row (linear interpolation)
function front(c, level = 0.5)
    i = findlast(>=(level), c)
    (i === nothing || i == nx) && return NaN
    return X[i] + (c[i] - level) / (c[i] - c[i + 1]) * dx
end
for (n, r) in enumerate(runs)
    @printf("Δρ = %2d: c = 0.5 front at top %.1f m, at bottom %.1f m\n", Δρ_ls[n], front(r.Cmat[1, :]),
        front(r.Cmat[ny, :]))
end
tol = 1.0e-12 # [m/s] Uy below this is round-off (Δρ = 0 gives |Uy| ~ 1e-17)

## plots
set_theme!(theme_latexfonts())
update_theme!(fontsize = 20)

skip = 4 # plot every 4th velocity vector
iq, jq = 2:skip:ny, 2:skip:nx
xq = vec([X[j] for i in iq, j in jq])
yq = vec([Ycoord[i] for i in iq, j in jq])
ascale = 3.5 / maximum(maxU) # same arrow scale in every panel [m per m/s]

function plume_axis(pos, title)
    return Axis(
        pos,
        title = title,
        xlabel = "x (m)",
        ylabel = "y (m), positive downward",
        yreversed = true, # y points down, as in the problem schematic
        aspect = DataAspect(),
        limits = (0, Lx, 0, Ly),
    )
end
function velocity_arrows!(ax, Ux, Uy; color = :white)
    arrows2d!(
        ax,
        xq,
        yq,
        vec(Ux[iq, jq]),
        vec(Uy[iq, jq]),
        lengthscale = ascale,
        color = color,
        shaftwidth = 1.5,
        tipwidth = 7,
        tiplength = 7,
    )
end

# Problem 1: concentration + Darcy velocity at day 7
fig = Figure(size = (900, 1250))
for (n, r) in enumerate(runs)
    ax = plume_axis(fig[n, 1], L"\Delta\rho = %$(Δρ_ls[n])\ \mathrm{kg/m^3},\ t = 7\ \mathrm{days}")
    heatmap!(ax, X, Ycoord, r.Cmat', colorrange = (0, 1), colormap = :viridis)
    contour!(ax, X, Ycoord, r.Cmat', levels = [0.5], color = :red, linewidth = 1.5)
    velocity_arrows!(ax, r.Ux, r.Uy)
end
Colorbar(fig[1:3, 2], limits = (0, 1), colormap = :viridis, label = "c (-)")
Label(
    fig[4, 1:2],
    @sprintf("red line: c = 0.5;  arrows: Darcy velocity (1 m of arrow = %.1e m/s)", 1 / ascale),
    fontsize = 15,
)
display(fig)
save("figures/hw6-p1-concentration-velocity.png", fig)

# Problem 3: upward-flow region for Δρ = 30
r30 = runs[3]
up = ifelse.(r30.Uy .< -tol, 1.0, NaN) # NaN cells are left transparent
fig3 = Figure(size = (900, 520))
ax3 = plume_axis(
    fig3[1, 1],
    L"\Delta\rho = 30\ \mathrm{kg/m^3},\ t = 7\ \mathrm{days}:\ U_y < 0\ \mathrm{(upward)\ shaded\ white}",
)
heatmap!(ax3, X, Ycoord, r30.Cmat', colorrange = (0, 1), colormap = :viridis)
heatmap!(ax3, X, Ycoord, up', colormap = [(:white, 0.5), (:white, 0.5)], nan_color = :transparent)
contour!(ax3, X, Ycoord, r30.Uy', levels = [0], color = :white, linewidth = 2.5)
velocity_arrows!(ax3, r30.Ux, r30.Uy, color = :black)
Colorbar(fig3[1, 2], limits = (0, 1), colormap = :viridis, label = "c (-)")
display(fig3)
save("figures/hw6-p3-upward.png", fig3)
