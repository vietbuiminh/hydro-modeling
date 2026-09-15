using LinearAlgebra
using CairoMakie

# physical properties
ρ = 1000 # g/cm^3
g = 9.81 # m/s^2
μ = 1.0e-3 # Pa*s

k = 1e-12 # m^2

# model dimension
Lx = 1
nx = 4
# nx = 1000 # uncomment to run the code for Problem 4
dx = Lx / nx

# boundaries conditions
u_in = 0.1e-2 / 60 # [m/s] specified Darcy flux at the left boundary (0.1 cm/min)
hR = 1.0 # [m] constant hydraulic head specified at the right boundary

# transmissibility
K = ρ*g*k/μ # m/s
Kmat = fill(K, nx, 1) # m/s

TL = fill(NaN, nx) # TL[1] unused: left boundary is flux-specified, not head-specified
TL[2:nx] = 2/dx * (1 ./ Kmat[1:(nx - 1)] + 1 ./ Kmat[2:nx]) .^ -1 # m/s

TR = fill(NaN, nx)
TR[nx] = 2/dx*Kmat[nx] # m/s
TR[1:(nx - 1)] = 2/dx * (1 ./ Kmat[1:(nx - 1)] + 1 ./ Kmat[2:nx]) .^ -1 # m/s

TC = fill(NaN, nx)
TC[1] = TR[1] # no left-transmissibility term: flux BC enters through b, not TC
TC[2:nx] = TL[2:nx] + TR[2:nx] # m/s

Tmat = Matrix(Tridiagonal(-TL[2:nx], TC, -TR[1:(nx - 1)])) # m/s

# RHS vector:b
b = zeros(nx)
b[1] = u_in # specified flux (Neumann)
b[nx] = TR[nx]*hR # specified head (Dirichlet)

# estimate h = T/b
h = Tmat \ b
xcell = range(dx/2, stop = Lx - dx/2, length = nx)

# analytical solution from Problem 1, h(x) = 1 + (u_in/K)(1-x)
h_analytical(x) = 1 .+ (u_in/K) .* (1 .- x)

# plot
set_theme!(theme_latexfonts())
update_theme!(
    fontsize = 20,
    figure_padding = (10, 50, 10, 10),
)

fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel = "x (m)",
    ylabel = "h (m)",
    title = "1D Discretization N = $nx",
)

x = range(0, stop = Lx, length = 100)

scatter!(ax, collect(xcell), h, markersize = 20, label = "h discretized")
lines!(ax, x, h_analytical.(x), label = "h analytical", color = :black)

axislegend(ax, position = :rt)
display(fig)
save("figures/1d-discretize-N$nx.png", fig)
