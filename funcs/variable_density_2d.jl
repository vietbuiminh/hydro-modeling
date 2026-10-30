using LinearAlgebra
using SparseArrays

### Variable-density flow coupled with advection-dispersion (Boussinesq)
# u = -K (∇h_f - β C e_y), β = Δρ/ρ, with +y pointing DOWN (y = 0 on the top boundary).
# Flow BCs: Neumann inflow uin [m/s] on the left, Dirichlet freshwater head hR on the
# right, no-flow top and bottom. Transport BCs: c = c_in(t) on the inflowing left
# boundary, zero gradient at the right outlet, no advective/dispersive flux top/bottom.
# Kmat is ny*nx hydraulic conductivity [m/s]; tsave are the times [s] to store.

# Transmissibility matrix: it does not depend on C, so it is factorized only once.
# The buoyancy term only enters the right-hand side.
function vd_flow_operator(Kmat, dx, dy)
    ny, nx = size(Kmat)
    nc = nx * ny

    # face conductivities, harmonic mean between neighbouring cells
    Kfx = 2 ./ (1 ./ Kmat[:, 1:(nx - 1)] .+ 1 ./ Kmat[:, 2:nx])
    Kfy = 2 ./ (1 ./ Kmat[1:(ny - 1), :] .+ 1 ./ Kmat[2:ny, :])

    TLmat = [zeros(ny, 1) (dy/dx) .* Kfx] # left boundary: Neumann, enters b
    TRmat = [(dy/dx) .* Kfx (2*dy/dx) .* Kmat[:, nx]] # right boundary: hR half a cell away
    TUmat = [zeros(1, nx); (dx/dy) .* Kfy] # top: no-flow
    TDmat = [(dx/dy) .* Kfy; zeros(1, nx)] # bottom: no-flow

    TL, TR, TU, TD = vec(TLmat), vec(TRmat), vec(TUmat), vec(TDmat)
    TC = TL .+ TR .+ TU .+ TD
    Tmat = spdiagm(
        0 => TC,
        -ny => -TL[(ny + 1):nc],
        -1 => -TU[2:nc],
        1 => -TD[1:(nc - 1)],
        ny => -TR[1:(nc - ny)],
    )
    return (; Tfac = lu(Tmat), Kfx, Kfy, TRmat)
end

# Head and face Darcy velocities for the current concentration field
function vd_flow(op, Cmat, dx, dy, uin, hR, β)
    (; Tfac, Kfx, Kfy, TRmat) = op
    ny, nx = size(Cmat)

    # buoyancy flux on the internal y faces (positive = downward), face C = arithmetic mean
    Cfy = 0.5 .* (Cmat[1:(ny - 1), :] .+ Cmat[2:ny, :])
    qb = Kfy .* β .* Cfy
    qbU = [zeros(1, nx); qb] # through the top face of each cell
    qbD = [qb; zeros(1, nx)] # through the bottom face

    # mass balance: T h = (buoyancy flux in - out) + boundary terms
    bmat = (qbU .- qbD) .* dx
    bmat[:, 1] .+= uin * dy # left boundary: specified inflow
    bmat[:, nx] .+= TRmat[:, nx] .* hR # right boundary: specified head
    hmat = reshape(Tfac \ vec(bmat), ny, nx)

    dh_x = hmat[:, 1:(nx - 1)] .- hmat[:, 2:nx]
    dh_y = hmat[1:(ny - 1), :] .- hmat[2:ny, :]
    ux = Kfx .* dh_x ./ dx # internal x faces
    uy = Kfy .* dh_y ./ dy .+ qb # internal y faces: -K (dh/dy - βC)

    UL = [fill(uin, ny, 1) ux]
    UR = [ux TRmat[:, nx] ./ dy .* (hmat[:, nx] .- hR)]
    UU = [zeros(1, nx); uy]
    UD = [uy; zeros(1, nx)]
    return (; hmat, UL, UR, UU, UD)
end

function solve_vd(Kmat, dx, dy, uin, hR, β, ϕ, D_diff, α, c_in, tsave)
    ny, nx = size(Kmat)
    op = vd_flow_operator(Kmat, dx, dy)

    Cmat = zeros(ny, nx)
    t = 0.0
    out = []
    for ts in tsave
        while t < ts
            (; UL, UR, UU, UD) = vd_flow(op, Cmat, dx, dy, uin, hR, β)

            # cell-centered pore velocity and dispersion coefficient (isotropic, as in class)
            Ux = 0.5 .* (UL .+ UR)
            Uy = 0.5 .* (UU .+ UD)
            V_mag = sqrt.(Ux .^ 2 .+ Uy .^ 2) ./ ϕ
            Dmat = D_diff .+ α .* V_mag

            Dx_har = 2 ./ (1 ./ Dmat[:, 1:(nx - 1)] .+ 1 ./ Dmat[:, 2:nx])
            Dy_har = 2 ./ (1 ./ Dmat[1:(ny - 1), :] .+ 1 ./ Dmat[2:ny, :])
            DL = [2 .* Dmat[:, 1] Dx_har] # inlet is half a cell away
            DR = [Dx_har zeros(ny, 1)] # zero-gradient outlet
            DU = [zeros(1, nx); Dy_har] # top: no flux
            DD = [Dy_har; zeros(1, nx)] # bottom: no flux

            # the velocity changes with C, so dt is recomputed every step
            dt_adv = 0.1 * min(dx, dy) / maximum(V_mag)
            dt_disp = 0.1 / (maximum(Dmat) * (1 / dx^2 + 1 / dy^2))
            dt = min(dt_adv, dt_disp, ts - t)

            CL = [fill(c_in(t), ny, 1) Cmat[:, 1:(nx - 1)]]
            CR = [Cmat[:, 2:nx] Cmat[:, nx]]
            CU = [zeros(1, nx); Cmat[1:(ny - 1), :]]
            CD = [Cmat[2:ny, :]; zeros(1, nx)]

            # advective flux, upwind
            Fadv_L = ifelse.(UL .>= 0, UL .* CL, UL .* Cmat)
            Fadv_R = ifelse.(UR .>= 0, UR .* Cmat, UR .* CR)
            Fadv_U = ifelse.(UU .>= 0, UU .* CU, UU .* Cmat)
            Fadv_D = ifelse.(UD .>= 0, UD .* Cmat, UD .* CD)
            netFadv = -Fadv_L .* dy .- Fadv_U .* dx .+ Fadv_D .* dx .+ Fadv_R .* dy

            # dispersive flux (positive = +x or +y)
            Fdis_L = ϕ .* DL .* (CL .- Cmat) ./ dx
            Fdis_R = ϕ .* DR .* (Cmat .- CR) ./ dx
            Fdis_U = ϕ .* DU .* (CU .- Cmat) ./ dy
            Fdis_D = ϕ .* DD .* (Cmat .- CD) ./ dy
            netFdis = -Fdis_L .* dy .- Fdis_U .* dx .+ Fdis_D .* dx .+ Fdis_R .* dy

            Cmat = Cmat .- dt / (ϕ * dx * dy) .* (netFadv .+ netFdis)
            t += dt
        end
        # flow consistent with the saved concentration field
        flow = vd_flow(op, Cmat, dx, dy, uin, hR, β)
        Ux = 0.5 .* (flow.UL .+ flow.UR)
        Uy = 0.5 .* (flow.UU .+ flow.UD)
        push!(out, (; t, Cmat = copy(Cmat), flow..., Ux, Uy))
    end
    return out
end
