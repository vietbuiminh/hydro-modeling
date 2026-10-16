using DelimitedFiles
using Statistics
using Printf
using CairoMakie

include("funcs/flow_transport_2d.jl")

# physical properties
ρ = 1000 # water density [kg/m^3]
μ = 1.0e-3 # dynamic viscosity [Pa s]
g = 9.81 # gravitational acceleration [m/s^2]
ϕ = 0.3 # porosity [-]

# model domain
Lx = 100 # domain length in x [m]
Ly = 50 # domain length in y [m]
nx = 100 # number of cells in x 
ny = 50 # number of cells in y
dx = Lx / nx # = 1 m
dy = Ly / ny # = 1 m

# boundary conditions
hL = 10 # left boundary head [m] (Dirichlet)
hR = 9 # right boundary head [m] (Dirichlet)
c_in = 0.0 # inflowing water
day = 24 * 60 * 60 # [s]

## get permeability fields
kmat = readdlm("data/permeability_field.csv", ',', Float64) # [m^2]
μ_log10k = -10 # mean of log10(k)
σ_log10k = 1 # std log10(k)
# field
Y = (log10.(kmat) .- μ_log10k) ./ σ_log10k

σ_ls = [0.1, 0.5, 1, 2] # heterogeneity levels
log10k(σ) = μ_log10k .+ σ .* Y

X = range(dx / 2, step = dx, length = nx) # cell-center coordinates [m]
Ycoord = range(dy / 2, step = dy, length = ny)

## initial plume: c = 1 for 0 <= x <= 10 m
C0 = zeros(ny, nx)
C0[:, X .<= 10] .= 1
M0 = ϕ * sum(C0) * dx * dy # initial solute mass (per unit thickness)

# 20,000 particles on a uniform 100 x 200 lattice over 0-10 m x 0-50 m
# (0.1 m x 0.25 m spacing, i.e. 40 particles per cell, 2,000 per 1-m column)
xp_ls = range(0.05, step = 0.1, length = 100)
yp_ls = range(0.125, step = 0.25, length = 200)
xp0 = vec([x for y in yp_ls, x in xp_ls])
yp0 = vec([y for y in yp_ls, x in xp_ls])
np = length(xp0)
np_per_m = np / 10 # = 2,000, normalizes the particle counts per 1-m interval

tplot = 5:5:30 # snapshot times [days]

# BTC time bins: 1 day wide for both methods. The run length grows with σ
# because low-k zones hold solute for a long time. For σ = 2 the particle
# arrivals reach ~3e7 days, far beyond what an explicit FVM can reach
# (dt ~ 2e-3 day), so that run is stopped at 20,000 days (~55 yr).
dtbin = 1 # [days]
tend = Dict(0.1 => 60, 0.5 => 250, 1 => 5000, 2 => 20000) # [days]

## all four heterogeneity levels
runs = map(σ_ls) do σ
    Kmat = ρ*g/μ .* 10 .^ log10k(σ) # hydraulic conductivity [m/s]
    flow = solve_flow(Kmat, dx, dy, hL, hR)

    tedges = collect(0:dtbin:tend[σ]) .* day
    fvm = solve_advection(flow, dx, dy, ϕ, C0, c_in, tedges, tplot .* day)
    ptcl = track_particles(flow, dx, dy, ϕ, xp0, yp0, tplot .* day)

    # P4-2: transverse-averaged concentration c̄(x,t) = (1/Ly) ∫ c dy
    cbar = [vec(sum(Cmat, dims = 1)) .* dy ./ Ly for Cmat in fvm.Cmat_ls]
    # particle counts in the same 1-m intervals, normalized by 2,000
    xbins = 0:dx:Lx
    pbar = map(eachcol(ptcl.xsnap)) do xs
        xin = filter(!isnan, xs)
        [count(x -> xbins[i] <= x < xbins[i + 1], xin) for i in 1:nx] ./ np_per_m
    end

    # P4-3: arrival-time PDFs [1/day] on the same bins
    tcen = (tedges[1:(end - 1)] .+ tedges[2:end]) ./ 2 ./ day
    f_fvm = fvm.Mout ./ sum(fvm.Mout) ./ dtbin # normalized by mass that exited
    tarr = ptcl.tarrive ./ day
    arrived = filter(<(tend[σ]), tarr)
    Ni = zeros(length(tcen)) # particle counts per time bin
    for t in arrived
        Ni[floor(Int, t / dtbin) + 1] += 1
    end
    f_ptcl = Ni ./ length(arrived) ./ dtbin # histogram with unit area

    (; σ, Kmat, flow, fvm, ptcl, cbar, pbar, tcen, f_fvm, f_ptcl,
        mass_out = sum(fvm.Mout) / M0, ptcl_out = length(arrived) / np,
        tmean_fvm = sum(tcen .* f_fvm) * dtbin, tmean_ptcl = mean(arrived),
        tmax_ptcl = maximum(tarr))
end

## report
for r in runs
    @printf("\nσ_log10k = %.1f (run to %d days)\n", r.σ, tend[r.σ])
    @printf("  max pore velocity = %.3f m/day, dt = %.4f day\n",
        maximum(r.fvm.V_mag) * day, r.fvm.dt_max / day)
    @printf("  FVM mass exited: %.4f, particles exited: %.4f\n",
        r.mass_out, r.ptcl_out)
    @printf("  mean arrival time (within run): FVM %.2f d, particles %.2f d\n",
        r.tmean_fvm, r.tmean_ptcl)
    @printf("  last particle arrival: %.3e days\n", r.tmax_ptcl)
end

## plots
set_theme!(theme_latexfonts())
update_theme!(fontsize = 18)

fieldaxis(pos, title) = Axis(
    pos,
    title = title,
    xlabel = "x (m)",
    ylabel = "y (m)",
    yreversed = true, # y points down, as in the problem schematic
    limits = (0, Lx, 0, Ly),
)

function rowlabels!(fig)
    for (i, r) in enumerate(runs)
        Label(
            fig[i, 1],
            L"\sigma_{\log_{10}k} = %$(r.σ)",
            rotation = π / 2,
            tellheight = false,
        )
    end
end

function hidegrid!(ax, i, j)
    i < length(runs) && hidexdecorations!(ax, grid = false) # x labels on bottom row
    j > 1 && hideydecorations!(ax, grid = false) # y labels on left column
end

# P4-1: FVM concentration fields, same color scale
fig = Figure(size = (1900, 800))
for (i, r) in enumerate(runs), (j, td) in enumerate(tplot)
    ax = fieldaxis(fig[i, j + 1], i == 1 ? "t = $td days" : "")
    heatmap!(
        ax,
        X,
        Ycoord,
        r.fvm.Cmat_ls[j]',
        colorrange = (0, 1),
        colormap = :viridis,
    )
    hidegrid!(ax, i, j)
end
rowlabels!(fig)
Colorbar(
    fig[:, length(tplot) + 2],
    limits = (0, 1),
    colormap = :viridis,
    label = "c (-)",
)
colgap!(fig.layout, 25) # keep neighbouring x tick labels ("100" and "0") apart
save("figures/ex-p4-concentration.png", fig)

# P4-1: particle positions at the same times, over the log10(k) field
krange = extrema(log10k(maximum(σ_ls)))
figp = Figure(size = (1900, 800))
for (i, r) in enumerate(runs), (j, td) in enumerate(tplot)
    ax = fieldaxis(figp[i, j + 1], i == 1 ? "t = $td days" : "")
    heatmap!(ax, X, Ycoord, log10k(r.σ)', colorrange = krange, colormap = :grays)
    scatter!(
        ax,
        r.ptcl.xsnap[:, j],
        r.ptcl.ysnap[:, j],
        markersize = 1,
        color = :red,
    )
    hidegrid!(ax, i, j)
end
rowlabels!(figp)
Colorbar(
    figp[:, length(tplot) + 2],
    limits = krange,
    colormap = :grays,
    label = L"\log_{10}k\ (m^2)",
)
colgap!(figp.layout, 25)
save("figures/ex-p4-particles.png", figp)

# P4-2: c̄(x,t) from FVM (lines) vs. normalized particle counts (steps)
colors = Makie.wong_colors()
figc = Figure(size = (1900, 500))
for (j, r) in enumerate(runs)
    ax = Axis(
        figc[1, j],
        title = L"\sigma_{\log_{10}k} = %$(r.σ)",
        xlabel = "x (m)",
        ylabel = L"\bar{c}(x,t)\ (-)",
        limits = (0, Lx, 0, nothing),
    )
    for (it, td) in enumerate(tplot)
        stairs!(ax, X, r.pbar[it], step = :center, color = (colors[it], 0.5))
        lines!(ax, X, r.cbar[it], color = colors[it], linewidth = 2.5,
            label = "t = $td d")
    end
    j > 1 && hideydecorations!(ax, grid = false, ticks = false, ticklabels = false)
end
Legend(
    figc[1, length(runs) + 1],
    [
        [LineElement(color = colors[it], linewidth = 2.5) for it in eachindex(tplot)];
        LineElement(color = :black, linewidth = 2.5);
        LineElement(color = (:black, 0.5))
    ],
    [["t = $td d" for td in tplot]; "FVM"; "particles"],
)
save("figures/ex-p4-cbar.png", figc)

# P4-3: outlet breakthrough PDFs, linear axes (top) and log-log for the tails
tshow = Dict(0.1 => 60, 0.5 => 200, 1 => 1000, 2 => 1000) # linear x-limits [days]
figb = Figure(size = (1900, 850))
for (j, r) in enumerate(runs)
    ax1 = Axis(
        figb[1, j],
        title = L"\sigma_{\log_{10}k} = %$(r.σ)",
        xlabel = "t (days)",
        ylabel = "pdf (1/day)",
        limits = (0, tshow[r.σ], 0, nothing),
    )
    stairs!(ax1, r.tcen, r.f_ptcl, step = :center, color = colors[2],
        label = "particles")
    lines!(ax1, r.tcen, r.f_fvm, color = colors[1], linewidth = 2.5, label = "FVM")

    ax2 = Axis(
        figb[2, j],
        xlabel = "t (days)",
        ylabel = "pdf (1/day)",
        xscale = log10,
        yscale = log10,
        limits = (1, tend[r.σ], 1e-7, 1),
    )
    pos = r.f_ptcl .> 0 # log axes: skip empty bins
    scatter!(ax2, r.tcen[pos], r.f_ptcl[pos], color = colors[2], markersize = 4)
    keep = r.f_fvm .> 1e-7
    lines!(ax2, r.tcen[keep], r.f_fvm[keep], color = colors[1], linewidth = 2.5)
    j == 1 && axislegend(ax1, position = :rt)
end
save("figures/ex-p4-btc.png", figb)
