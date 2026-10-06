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
hL = 10 # left boundary head [m] (Dirichlet)
hR = 9 # right boundary head [m] (Dirichlet)
day = 24 * 60 * 60 # [s]
c_in(t) = t < 10 * day ? 1.0 : 0.0 # inject c = 1 for the first 10 days

## flow
Kmat = (ρ*g*k/μ) .* ones(ny, nx) # hydraulic conductivity [m/s]
flow = solve_flow(Kmat, dx, dy, hL, hR)

## transport, saved every day so the plume can be tracked through time
tdays = 1:30
res = solve_transport(flow, dx, dy, ϕ, D_diff, α, c_in, tdays .* day)

X = range(dx / 2, step = dx, length = nx) # cell-center coordinates [m]
Ycoord = range(dy / 2, step = dy, length = ny)

# transverse-averaged concentration c̄(x,t) = (1/Ly) ∫ c dy
cbar = [vec(sum(Cmat, dims = 1)) .* dy ./ Ly for Cmat in res.Cmat_ls]

## expected mean pore velocity from the total inflow at the left boundary
ux_mean = sum(flow.UL[:, 1]) * dy / Ly # [m/s]
vx_mean = ux_mean / ϕ # [m/s]
vx_day = vx_mean * day # [m/day]

# plume centroid and position of the c̄ = 0.5 leading front
centroid(c) = sum(X .* c) / sum(c)
function front(c, level = 0.5)
    i = findlast(>=(level), c)
    (i === nothing || i == nx) && return NaN
    return X[i] + (c[i] - level) / (c[i] - c[i + 1]) * dx # linear interpolation
end
xbar = centroid.(cbar)
xfront = front.(cbar)

# expected positions for a 10-day slug moving at vx: front at vx*t, centroid at
# the middle of the injected slug
xfront_exp = vx_day .* tdays
xbar_exp = [t <= 10 ? vx_day * t / 2 : vx_day * (t - 5) for t in tdays]

# plume speeds: front during days 5-20, centroid once injection has stopped
# and before the plume reaches the outlet (days 12-22)
v_front = (xfront[20] - xfront[5]) / (20 - 5) # [m/day]
v_centroid = (xbar[22] - xbar[12]) / (22 - 12) # [m/day]

# solute mass in the domain vs. mass that entered by advection
mass = [ϕ * sum(Cmat) * dx * dy for Cmat in res.Cmat_ls] # [kg/m per unit c]
mass_in = ux_mean * Ly * 10 * day # advective inflow over 10 days

## report
Cmin = minimum(minimum.(res.Cmat_ls))
Cmax = maximum(maximum.(res.Cmat_ls))
transverse_dev =
    maximum(maximum(abs.(Cmat .- mean(Cmat, dims = 1))) for Cmat in res.Cmat_ls)

@printf("K = %.3e m/s\n", Kmat[1])
@printf("head: %.4f to %.4f m\n", extrema(flow.hmat)...)
@printf("mean Darcy velocity ux = %.4e m/s\n", ux_mean)
@printf("mean pore velocity vx = %.4e m/s = %.4f m/day\n", vx_mean, vx_day)
@printf("max |Uy| = %.3e m/s\n", maximum(abs.(0.5 .* (flow.UU .+ flow.UD))))
@printf("D_h = %.4e m^2/s\n", maximum(res.Dmat))
@printf("dt_adv = %.1f s, dt_disp = %.1f s\n", res.dt_adv, res.dt_disp)
@printf("c range over all saved times: [%.3e, %.6f]\n", Cmin, Cmax)
@printf("max transverse deviation |c - c̄|: %.3e\n", transverse_dev)
@printf("front speed (day 5-20) = %.4f m/day\n", v_front)
@printf("centroid speed (day 12-22) = %.4f m/day\n", v_centroid)
println("\n t[d]  x_front  vx*t   x̄      x̄_exp   max c̄   mass/mass_in")
for i in tdays
    @printf("%4d  %6.2f  %6.2f  %6.2f  %6.2f  %.3f  %.3f\n", i, xfront[i],
        xfront_exp[i], xbar[i], xbar_exp[i], maximum(cbar[i]), mass[i] / mass_in)
end

## plots
set_theme!(theme_latexfonts())
update_theme!(fontsize = 18)
tplot = 5:5:30

# 2-D concentration fields
fig = Figure(size = (1080, 600))
for (n, td) in enumerate(tplot)
    r, col = fldmod1(n, 3)
    ax = Axis(
        fig[r, col],
        title = "t = $td days",
        xlabel = "Lx (m)",
        ylabel = "Ly (m)",
        yreversed = true, # y points down, as in the problem schematic
        # aspect = DataAspect(),
    )
    heatmap!(ax, X, Ycoord, res.Cmat_ls[td]', colorrange = (0, 1), colormap = :viridis)
end
Colorbar(fig[:, 4], limits = (0, 1), colormap = :viridis, label = "c (-)")
save("figures/hw5-p1-concentration.png", fig)

# transverse-averaged concentration profiles
fig2 = Figure(size = (900, 500))
ax2 = Axis(
    fig2[1, 1],
    xlabel = "Lx (m)",
    ylabel = L"\bar{c}(x,t)\ (-)",
    limits = (0, Lx, 0, 1.05),
)
for td in tplot
    lines!(ax2, X, cbar[td], label = "t = $td d", linewidth = 2)
end
axislegend(ax2, position = :rt)
save("figures/hw5-p1-cbar.png", fig2)

# fig2_gif = Figure(size = (900, 500))
# ax2_gif = Axis(fig2_gif[1, 1], xlabel = "Lx (m)", ylabel = L"\bar{c}(x,t)\ (-)",
# limits = (0, Lx, 0, 1.05))
# record(fig2_gif, "figures/hw5-p1-cbar.gif", tplot; framerate = 2) do td
#     lines!(ax2_gif, X, cbar[td], label = "t = $td d", linewidth = 2)
# end

# plume position vs. time compared with the mean pore velocity
fig3 = Figure(size = (900, 500))
ax3 = Axis(fig3[1, 1], xlabel = "t (days)", ylabel = "position (m)")
lines!(ax3, tdays, xfront_exp, color = :gray, linestyle = :dash, label = L"\bar{v}_x t")
scatter!(ax3, tdays, xfront, label = L"front, $\bar{c} = 0.5$")
lines!(
    ax3,
    tdays,
    xbar_exp,
    color = :black,
    linestyle = :dot,
    label = "expected centroid",
)
scatter!(ax3, tdays, xbar, marker = :diamond, label = L"centroid $\bar{x}$")
axislegend(ax3, position = :lt)
save("figures/hw5-p1-plume-position.png", fig3)
