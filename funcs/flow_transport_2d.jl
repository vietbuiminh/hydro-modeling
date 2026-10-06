using LinearAlgebra

### Steady-state 2-D flow (TPFA) 
# Dirichlet heads hL, hR on the L/R boundaries and no-flow top/bottom.
# Kmat is ny*nx hydraulic conductivity [m/s]. 
# Returns the head field and the Darcy velocities [m/s] on
# the L/R/U/D faces.
function solve_flow(Kmat, dx, dy, hL, hR)
    ny, nx = size(Kmat)
    nc = nx * ny

    K_har_x = (1 ./ Kmat[:, 1:(nx - 1)] .+ 1 ./ Kmat[:, 2:nx]) .^ -
    1 # harmonic mean in x direction
    K_har_y = (1 ./ Kmat[1:(ny - 1), :] .+ 1 ./ Kmat[2:ny, :]) .^ -
    1 # harmonic mean in y direction

    TLmat = fill(NaN, ny, nx)
    TLmat[:, 1] .= (2*dy/dx) .* Kmat[:, 1] # left boundary (Dirichlet), half-
    cell distance
    TLmat[:, 2:nx] .= (2*dy/dx) .* K_har_x

    TRmat = fill(NaN, ny, nx)
    TRmat[:, nx] .= (2*dy/dx) .* Kmat[:, nx] # right boundary (Dirichlet), half-
    cell distance
    TRmat[:, 1:(nx - 1)] .= (2*dy/dx) .* K_har_x

    TUmat = fill(NaN, ny, nx)
    TUmat[1, :] .= 0.0 # top boundary: no-flow
    TUmat[2:ny, :] .= (2*dx/dy) .* K_har_y

    TDmat = fill(NaN, ny, nx)
    TDmat[ny, :] .= 0.0 # bottom boundary: no-flow
    TDmat[1:(ny - 1), :] .= (2*dx/dy) .* K_har_y

    TL = reshape(TLmat, nc)
    TR = reshape(TRmat, nc)
    TU = reshape(TUmat, nc)
    TD = reshape(TDmat, nc)
    TC = TL .+ TR .+ TU .+ TD

    Tmat =
        diagm(0 => TC) +
        diagm(-ny => -TL[(ny + 1):nc]) +
        diagm(-1 => -TU[2:nc]) +
        diagm(1 => -TD[1:(nc - 1)]) +
        diagm(ny => -TR[1:(nc - ny)])

    bmat = zeros(ny, nx)
    bmat[:, 1] .= TLmat[:, 1] .* hL # left boundary: specified head
    bmat[:, nx] .= TRmat[:, nx] .* hR # right boundary: specified head
    hmat = reshape(Tmat \ reshape(bmat, nc), ny, nx)

    # Darcy velocity on faces: transmissibility / face length * head difference
    dh_x = hmat[:, 1:(nx - 1)] .- hmat[:, 2:nx]
    dh_y = hmat[1:(ny - 1), :] .- hmat[2:ny, :]

    UL = fill(NaN, ny, nx)
    UL[:, 1] .= TLmat[:, 1] ./ dy .* (hL .- hmat[:, 1])
    UL[:, 2:nx] .= TLmat[:, 2:nx] ./ dy .* dh_x

    UR = fill(NaN, ny, nx)
    UR[:, nx] .= TRmat[:, nx] ./ dy .* (hmat[:, nx] .- hR)
    UR[:, 1:(nx - 1)] .= TRmat[:, 1:(nx - 1)] ./ dy .* dh_x

    UU = fill(NaN, ny, nx)
    UU[1, :] .= 0.0
    UU[2:ny, :] .= TUmat[2:ny, :] ./ dx .* dh_y

    UD = fill(NaN, ny, nx)
    UD[ny, :] .= 0.0
    UD[1:(ny - 1), :] .= TDmat[1:(ny - 1), :] ./ dx .* dh_y

    return (; hmat, UL, UR, UU, UD)
end

# Explicit (forward Euler) advection-dispersion with 1st-order upwinding.
# tsave: the times [s] at which to store the concentration field; the time
# step is shortened so the simulation lands exactly on each of them.
function solve_transport(flow, dx, dy, ϕ, D_diff, α, c_in, tsave)
    (; UL, UR, UU, UD) = flow
    ny, nx = size(UL)

    # cell-centered pore velocity and hydrodynamic dispersion coefficient
    Ux = 0.5 .* (UL .+ UR)
    Uy = 0.5 .* (UU .+ UD)
    V_mag = sqrt.(Ux .^ 2 .+ Uy .^ 2) ./ ϕ
    Dmat = D_diff .+ α .* V_mag

    # D_h on faces: harmonic mean inside, cell value over the half cell at the inlet
    Dx_har = 2 ./ (1 ./ Dmat[:, 1:(nx - 1)] .+ 1 ./ Dmat[:, 2:nx])
    Dy_har = 2 ./ (1 ./ Dmat[1:(ny - 1), :] .+ 1 ./ Dmat[2:ny, :])
    DL = [2 .* Dmat[:, 1] Dx_har] # 2D/dx == D/(dx/2): inlet is half a cell away
    DR = [Dx_har zeros(ny, 1)] # zero-gradient outlet: no dispersive flux
    DU = [zeros(1, nx); Dy_har] # top: no flux
    DD = [Dy_har; zeros(1, nx)] # bottom: no flux

    # dt is fixed because the velocity field is steady
    dt_adv = 0.1 * min(dx, dy) / maximum(V_mag)
    dt_disp = 0.1 / (maximum(Dmat) * (1 / dx^2 + 1 / dy^2))
    dt_max = min(dt_adv, dt_disp)

    Cmat = zeros(ny, nx)
    t = 0.0
    Cmat_ls = Matrix{Float64}[] # storage to return the field later
    for ts in tsave
        while t < ts
            dt = min(dt_max, ts - t)

            # cell-centered concentration
            CL = [fill(c_in(t), ny, 1) Cmat[:, 1:(nx - 1)]]
            CR = [Cmat[:, 2:nx] Cmat[:, nx]]
            CU = [zeros(1, nx); Cmat[1:(ny - 1), :]]
            CD = [Cmat[2:ny, :]; zeros(1, nx)]

            # advective flux, upwind
            Fadv_L = ifelse.(UL .>= 0, UL .* CL, UL .* Cmat)
            Fadv_R = ifelse.(UR .>= 0, UR .* Cmat, UR .* CR)
            Fadv_U = ifelse.(UU .>= 0, UU .* CU, UU .* Cmat)
            Fadv_D = ifelse.(UD .>= 0, UD .* Cmat, UD .* CD)

            # net flux advection
            netFadv = -Fadv_L .* dy .- Fadv_U .* dx .+ Fadv_D .* dx .+ Fadv_R .* dy

            # dispersive flux, cell-centered differences (positive = +x or +y)
            Fdis_L = ϕ .* DL .* (CL .- Cmat) ./ dx
            Fdis_R = ϕ .* DR .* (Cmat .- CR) ./ dx
            Fdis_U = ϕ .* DU .* (CU .- Cmat) ./ dy
            Fdis_D = ϕ .* DD .* (Cmat .- CD) ./ dy

            # net flux dispersion
            netFdis = -Fdis_L .* dy .- Fdis_U .* dx .+ Fdis_D .* dx .+ Fdis_R .* dy

            Cmat = Cmat .- dt / (ϕ * dx * dy) .* (netFadv + netFdis)
            t += dt
        end
        push!(Cmat_ls, copy(Cmat))
    end
    return (; Cmat_ls, V_mag, Dmat, dt_adv, dt_disp)
end
