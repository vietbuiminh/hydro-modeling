using DelimitedFiles
using Statistics
using Printf
using CairoMakie

include("funcs/flow_transport_2d.jl")

# HW5 Problem 2: advection-dispersion in heterogeneous aquifers

# physical properties
ρ = 1000 # water density [kg/m^3]
μ = 1.0e-3 # dynamic viscosity [Pa s]
g = 9.81 # gravitational acceleration [m/s^2]

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

## permeability fields
kmat = readdlm("data/permeability_field.csv", ',', Float64) # [m^2]
μ_log10k = -10 # given mean of log10(k)
σ_log10k = 1 # given standard deviation of log10(k)
# standardized field, same pattern for every σ
Y = (log10.(kmat) .- μ_log10k) ./ σ_log10k

σ_ls = [0.1, 0.5, 1, 2] # heterogeneity levels a.k.a. standard deviations of log10(k)
log10k(σ) = μ_log10k .+ σ .* Y

X = range(dx / 2, step = dx, length = nx) # cell-center coordinates [m]
Ycoord = range(dy / 2, step = dy, length = ny)
tdays = 1:30 # saved every day so x̄ and σx can be tracked through time
tplot = 5:5:30

centroid(c) = sum(X .* c) / sum(c)
variance(c) = sqrt(sum((X .- centroid(c)) .^ 2 .* c) / sum(c)) # σx

# σx at a given centroid location, by linear interpolation between saved times
function σx_at(xb, sx, x0)
    i = findfirst(>=(x0), xb)
    (i === nothing || i == 1) && return NaN
    w = (x0 - xb[i - 1]) / (xb[i] - xb[i - 1])
    return (1 - w) * sx[i - 1] + w * sx[i], (i - 1) + w # σx, time [days]
end

## all four heterogeneity levels
runs = map(σ_ls) do σ
    Kmat = ρ*g/μ .* 10 .^ log10k(σ) # hydraulic conductivity [m/s]
    flow = solve_flow(Kmat, dx, dy, hL, hR)
    res = solve_transport(flow, dx, dy, ϕ, D_diff, α, c_in, tdays .* day)

    # transverse-averaged concentration c̄(x,t) = (1/Ly) ∫ c dy
    cbar = [vec(sum(Cmat, dims = 1)) .* dy ./ Ly for Cmat in res.Cmat_ls]
    ux_mean = sum(flow.UL[:, 1]) * dy / Ly # mean Darcy velocity from the inflow [m/s]
    mass = [ϕ * sum(Cmat) * dx * dy for Cmat in res.Cmat_ls]
    mass_in = ux_mean * Ly * 10 * day # advective inflow over 10 days

    (; σ, Kmat, flow, res, cbar, ux_mean,
        vx_day = ux_mean / ϕ * day, # mean pore velocity [m/day]
        # flux-weighted inlet pore velocity: solute enters where the inflow is fastest
        vfw_day = sum(flow.UL[:, 1] .^ 2) / sum(flow.UL[:, 1]) / ϕ * day,
        xbar = centroid.(cbar), σx = variance.(cbar),
        mass_frac = mass ./ mass_in,
        cmin = minimum(minimum.(res.Cmat_ls)), cmax = maximum(maximum.(res.Cmat_ls)))
end

## report
# for r in runs
#     @printf("\nσ_log10k = %.1f\n", r.σ)
#     @printf("  K range: %.3e to %.3e m/s\n", extrema(r.Kmat)...)
#     @printf("  mean pore velocity vx = %.4f m/day\n", r.vx_day)
#     @printf("  max |v| = %.4f m/day\n", maximum(r.res.V_mag) * day)
#     @printf("  dt_adv = %.1f s, dt_disp = %.1f s\n", r.res.dt_adv, r.res.dt_disp)
#     @printf("  c range: [%.3e, %.6f]\n", r.cmin, r.cmax)
#     @printf("  flux-weighted inlet pore velocity = %.4f m/day\n", r.vfw_day)
#     @printf(
#         "  centroid speed during injection, 2 x̄(5 d)/5 d = %.4f m/day\n",
#         2 * r.xbar[5] / 5
#     )
#     sx50, t50 = σx_at(r.xbar, r.σx, 50)
#     i = floor(Int, t50)
#     m50 = r.mass_frac[i] + (t50 - i) * (r.mass_frac[i + 1] - r.mass_frac[i])
#     @printf(
#         "  at x̄ = 50 m: σx = %.2f m, t = %.2f days, mass/mass_in = %.3f\n",
#         sx50,
#         t50,
#         m50
#     )
#     println("  t[d]   x̄      σx    mass/mass_in")
#     for td in tplot
#         @printf(
#             "  %3d  %6.2f  %6.2f  %.3f\n",
#             td,
#             r.xbar[td],
#             r.σx[td],
#             r.mass_frac[td]
#         )
#     end
# end

## plots
set_theme!(theme_latexfonts())
update_theme!(fontsize = 18)

fieldaxis(pos, title) = Axis(
    pos,
    title = title,
    xlabel = "x (m)",
    ylabel = "y (m)",
    yreversed = true, # y points down, as in the problem schematic
    # aspect = DataAspect(),
)

# 2D concentration fields, same color scale
fig = Figure(size = (1900, 800))
for (i, r) in enumerate(runs), (j, td) in enumerate(tplot)
    ax = fieldaxis(fig[i, j + 1], i == 1 ? "t = $td days" : "")
    heatmap!(
        ax,
        X,
        Ycoord,
        r.res.Cmat_ls[td]',
        colorrange = (0, 1),
        colormap = :viridis,
    )
    i < length(runs) && hidexdecorations!(ax, grid =
    false) # x labels only on the bottom row
    j > 1 && hideydecorations!(ax, grid = false) # y labels only on the left column
end
for (i, r) in enumerate(runs)
    Label(
        fig[i, 1],
        L"\sigma_{\log_{10}k} = %$(r.σ)",
        rotation = π / 2,
        tellheight = false,
    )
end
Colorbar(
    fig[:, length(tplot) + 2],
    limits = (0, 1),
    colormap = :viridis,
    label = "c (-)",
)
colgap!(fig.layout, 25) # keep neighbouring x tick labels ("100" and "0") apart
save("figures/hw5-p2-concentration.png", fig)

# permeability fields (top row, shared color scale) and transverse-averaged
# concentration profiles (bottom row), one column per σ
krange = extrema(log10k(maximum(σ_ls)))
figkc = Figure(size = (1900, 750))
for (j, r) in enumerate(runs)
    axk = fieldaxis(figkc[1, j], L"\sigma_{\log_{10}k} = %$(r.σ)")
    heatmap!(axk, X, Ycoord, log10k(r.σ)', colorrange = krange, colormap = :viridis)

    axc = Axis(
        figkc[2, j],
        xlabel = "x (m)",
        ylabel = L"\bar{c}(x,t)\ (-)",
        limits = (0, Lx, 0, 1.05),
    )
    for td in tplot
        lines!(axc, X, r.cbar[td], label = "t = $td d", linewidth = 2)
    end
    if j > 1
        hideydecorations!(axk, grid = false) # y labels only on the left column
        hideydecorations!(axc, grid = false)
    end
end
Colorbar(
    figkc[1, length(runs) + 1],
    limits = krange,
    colormap = :viridis,
    label = L"\log_{10}k\ (m^2)",
)
Legend(figkc[2, length(runs) + 1], figkc.content[2])
colgap!(figkc.layout, 25)
save("figures/hw5-p2-log10k-cbar.png", figkc)

# centroid through time and spreading vs. centroid
figm = Figure(size = (1200, 500))
ax1 = Axis(figm[1, 1], xlabel = "t (days)", ylabel = L"centroid $\bar{x}$ (m)")
ax2 = Axis(figm[1, 2], xlabel = L"centroid $\bar{x}$ (m)", ylabel = L"\sigma_x\ (m)")
for r in runs
    lab = L"\sigma_{\log_{10}k} = %$(r.σ)"
    scatterlines!(ax1, tdays, r.xbar, label = lab)
    scatterlines!(ax2, r.xbar, r.σx, label = lab)
end
vlines!(ax2, [50], color = :gray, linestyle = :dash)
axislegend(ax1, position = :lt)
axislegend(ax2, position = :lt)
save("figures/hw5-p2-moments.png", figm)
