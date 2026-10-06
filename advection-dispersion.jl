using LinearAlgebra
using DelimitedFiles
using Statistics
using CairoMakie

# physical properties
ρ = 1000 # water density [kg/m^3]
μ = 1.0e-3 # dynamic viscosity [Pa s]
g = 9.81 # gravitational acceleration [m/s^2]

ϕ = 0.3 # porosity [-]
D_diff = 1.0e-9 # molecular diffusion coeff. [m^2/s]
α = 1 # isotropic dispersivity [m]

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
σ = 0.5 # heterogeneity level
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

## pore velocity
V_mag = sqrt.(Ux .^ 2 .+ Uy .^ 2) ./ ϕ # field pore-velocity magnitude
V_max = maximum(V_mag)

## transport (advection + dispersion)
# the time loop lives in a `let` block so its variables are local; at the top
# level of a script, assigning globals inside a `while` loop is ambiguous
Cmat = let
    t = 0.0 # time [s]
    Cmat = zeros(ny, nx)
    C_in = 1.0 # left boundary concentration (Dirichlet)

    doPlot = false # should I plot?
    tPlot = (1:30) .* 24*60*60 # when to plot
    itPlot = 1 # counting plots

    # concentration figure: Observables let us update it in place, like drawnow
    Cobs = Observable(collect(Cmat'))
    titleobs = Observable("time: 0.00 days")
    figC = Figure(size = (800, 1000))
    axK = Axis(figC[1, 1], title = "log10(K)")
    heatmap!(axK, X, Ycoord, log10.(Kmat)')
    axH = Axis(figC[2, 1], title = "Hydraulic head")
    heatmap!(axH, X, Ycoord, hmat')
    add_velocity!(axH)
    axC = Axis(figC[3, 1], xlabel = "x [m]", ylabel = "y [m]", title = titleobs)
    heatmap!(axC, X, Ycoord, Cobs, colorrange = (0, 1))

    while t < 30*24*60*60
        if t > 5*24*60*60
            C_in = 0.0 # stop injecting after 5 days
        end

        ## neighboring cell centers' concentration
        CL = [fill(C_in, ny, 1) Cmat[:, 1:(nx - 1)]]
        CR = [Cmat[:, 2:nx] Cmat[:, nx]] # zero-
        gradient concentration condition [Neumann]
        CU = [zeros(1, nx); Cmat[1:(ny - 1), :]]
        CD = [Cmat[2:ny, :]; zeros(1, nx)]

        ## advective flux estimation (on cell faces)
        Fadv_L = fill(NaN, ny, nx)
        Fadv_U = fill(NaN, ny, nx)
        Fadv_D = fill(NaN, ny, nx)
        Fadv_R = fill(NaN, ny, nx)

        # upwinding
        # left
        ind = UL .>= 0
        Fadv_L[ind] .= UL[ind] .* CL[ind]
        Fadv_L[.!ind] .= UL[.!ind] .* Cmat[.!ind]

        # right
        ind = UR .>= 0
        Fadv_R[ind] .= UR[ind] .* Cmat[ind]
        Fadv_R[.!ind] .= UR[.!ind] .* CR[.!ind]

        # up
        ind = UU .>= 0
        Fadv_U[ind] .= UU[ind] .* CU[ind]
        Fadv_U[.!ind] .= UU[.!ind] .* Cmat[.!ind]

        # down
        ind = UD .>= 0
        Fadv_D[ind] .= UD[ind] .* Cmat[ind]
        Fadv_D[.!ind] .= UD[.!ind] .* CD[.!ind]

        # net flux
        netFadv = -Fadv_L .* dy .- Fadv_U .* dx .+ Fadv_D .* dx .+ Fadv_R .* dy

        ## dispersive flux estimation (on cell faces)
        # dispersion coeff on cell centers (D: D_diff + α * V_mag)
        Dmat = D_diff .+ α .* V_mag

        # dispersion coeff on internal faces (excluding domain boundary faces)
        Dx_har = 2 ./ (1 ./ Dmat[:, 1:(nx - 1)] .+ 1 ./ Dmat[:, 2:nx])
        Dy_har = 2 ./ (1 ./ Dmat[1:(ny - 1), :] .+ 1 ./ Dmat[2:ny, :])

        # dispersion coeff on cell faces
        DL = fill(NaN, ny, nx)
        DR = fill(NaN, ny, nx)
        DU = fill(NaN, ny, nx)
        DD = fill(NaN, ny, nx)

        DL[:, 1] .= 2 .* Dmat[:, 1]
        DL[:, 2:nx] .= Dx_har

        DR[:, nx] .= Dmat[:, nx]
        DR[:, 1:(nx - 1)] .= Dx_har

        DU[1, :] .= 0.0
        DU[2:ny, :] .= Dy_har

        DD[ny, :] .= 0.0
        DD[1:(ny - 1), :] .= Dy_har

        # flux
        Fdis_L = ϕ .* DL .* (CL .- Cmat) ./ dx
        Fdis_R = ϕ .* DR .* (Cmat .- CR) ./ dx
        Fdis_U = ϕ .* DU .* (CU .- Cmat) ./ dy
        Fdis_D = ϕ .* DD .* (Cmat .- CD) ./ dy

        # net flux
        netFdis = -Fdis_L .* dy .- Fdis_U .* dx .+ Fdis_D .* dx .+ Fdis_R .* dy

        ## determine dt
        dt_adv = 0.1 * dx / V_max
        dt_dis = 0.1 * dx^2 / maximum(Dmat)
        dt = min(dt_adv, dt_dis)

        if itPlot <= length(tPlot) && t + dt > tPlot[itPlot]
            dt = tPlot[itPlot] - t
            doPlot = true
            itPlot += 1
        end

        ## advance concentration field in time
        t += dt
        Cmat = Cmat .- dt / (ϕ * dx * dy) .* (netFadv .+ netFdis)

        if doPlot
            ## visualization
            Cobs[] = collect(Cmat')
            titleobs[] = "time: $(round(t / 24 / 60 / 60, digits = 2)) days"
            display(figC)
            doPlot = false
            sleep(0.1) # pause to let the figure update
        end
    end
    Cmat
end
