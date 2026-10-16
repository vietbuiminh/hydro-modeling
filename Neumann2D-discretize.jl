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
hL = 10 # left boundary head [m] (Dirichlet type)
hR = 9 # right boundary head [m] (Dirichlet type)

# read in permeability field from CSV file
kmat = readdlm("data/permeability_field.csv", ',', Float64) # [m^2]
Kmat = ρ*g*kmat/μ .* kmat
kmean = mean(log10.(reshape(kmat, nc))) # log10 needs the broadcast dot for a Vector
kstd = std(log10.(reshape(kmat, nc)))

Y = (log10.(kmat) .- kmean) ./ kstd # standardized spatial field

σ_ls = [0.1, 0.5, 1, 2] # heterogeneity levels to compare

set_theme!(theme_latexfonts())
update_theme!(
    fontsize = 18,
    Axis = (yreversed = true,),
)
X = range(dx / 2, step = dx, length = nx)
Ycoord = range(dy / 2, step = dy, length =
ny) # renamed: `Y` above is the standardized field, don't overwrite it

# this could be a function, but for the sake of the hw and deadline, I am leaving it as
# a loop
# TODO (in the future): make a function :)
for σ in σ_ls
    log10k = kmean .+ σ .* Y # rescaled log10(k) field
    Kmat = ρ*g/μ .* 10 .^ log10k # hydraulic conductivity [m/s]

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
    TUmat[2:ny, :] .= (2*dx/dy) .* K_har_y # transmissibility in y direction

    TDmat = fill(NaN, ny, nx)
    TDmat[ny, :] .= 0.0 # bottom boundary: no transmissibility
    TDmat[1:(ny - 1), :] .= (2*dx/dy) .* K_har_y # transmissibility in y direction

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

    # hydraulic head vector h
    h = Tmat \ b # solve for hydraulic head
    hmat = reshape(h, ny, nx) # reshape h vector to matrix form

    # Darcy velocity at cell faces
    dh_x = hmat[:, 1:(nx - 1)] .- hmat[:, 2:nx] # head difference in x direction
    dh_y = hmat[1:(ny - 1), :] .- hmat[2:ny, :] # head difference in y direction

    UL = fill(NaN, ny, nx)
    UL[:, 1] .= TLmat[:, 1] .* dh_x[:,
        1] ./ dx # left boundary (Dirichlet): nearest interior gradient
    UL[:, 2:nx] .= TLmat[:, 2:nx] .* dh_x ./ dx # Darcy velocity in x direction

    UR = fill(NaN, ny, nx)
    UR[:, nx] .= TRmat[:, nx] .* dh_x[:, nx - 1] ./ dx # right boundary (Dirichlet): nearest interior gradient
    UR[:, 1:(nx - 1)] .= TRmat[:, 1:(nx - 1)] .*
                         dh_x ./ dx # Darcy velocity in x direction

    UU = fill(NaN, ny, nx)
    UU[1, :] .= 0.0 # top boundary (Dirichlet): no-flow
    UU[2:ny, :] .= TUmat[2:ny, :] .* dh_y ./ dy # Darcy velocity in y direction

    UD = fill(NaN, ny, nx)
    UD[ny, :] .= 0.0 # bottom boundary (Dirichlet): no-flow
    UD[1:(ny - 1), :] .= TDmat[1:(ny - 1), :] .*
                         dh_y ./ dy # Darcy velocity in y direction

    Ux = 0.5 * (UL .+ UR) # cell-centered Darcy velocity, x direction
    Uy = 0.5 * (UU .+ UD) # cell-centered Darcy velocity, y direction

    # plot
    fig = Figure(size = (1150, 420))

    ax1 = Axis(
        fig[1, 1],
        xlabel = "Lx (m)", ylabel = "Ly (m)", title = "log10(k)")
    hm1 = heatmap!(ax1, X, Ycoord, log10k')
    Colorbar(fig[1, 2], hm1)

    ax2 = Axis(
        fig[1, 3],
        xlabel = "Lx (m)",
        ylabel = "Ly (m)",
        title = "Hydraulic Head & Darcy Velocity (σ = $σ)",
    )
    hm2 = heatmap!(ax2, X, Ycoord, hmat')
    Colorbar(fig[1, 4], hm2, label = "(m)")

    # overlay velocity vectors on the hydraulic head heatmap (subsample for clarity)
    xCells = 1:nx # cell-index coordinates, matching the heatmap's default axes
    yCells = 1:ny
    skip = 5

    # Darcy velocities are tiny relative to the ~100-cell axes, so the raw
    # vectors are invisible; lengthscale blows them up to a few cells long
    maxU = maximum(sqrt.(Ux .^ 2 .+ Uy .^ 2))
    arrows2d!(
        ax2,
        xCells[1:skip:end],
        yCells[1:skip:end],
        Ux[1:skip:end, 1:skip:end]',
        Uy[1:skip:end, 1:skip:end]',
        color = :black,
        lengthscale = 2 * skip / maxU,
    )

    display(fig)
    save("figures/hw4-p3-sigma$(σ).png", fig)
end
