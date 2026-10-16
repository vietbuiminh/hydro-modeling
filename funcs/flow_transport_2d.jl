using LinearAlgebra

### Steady-state 2-D flow (TPFA) 
# Dirichlet heads hL, hR on the L/R boundaries and no-flow top/bottom.
# Kmat is ny*nx hydraulic conductivity [m/s]. 
# Returns the head field and the Darcy velocities [m/s] on
# the L/R/U/D faces.
function solve_flow(Kmat, dx, dy, hL, hR)
    ny, nx = size(Kmat)
    nc = nx * ny

    # harmonic means in x and y directions
    K_har_x = 1 ./ (1 ./ Kmat[:, 1:(nx - 1)] .+ 1 ./ Kmat[:, 2:nx])
    K_har_y = 1 ./ (1 ./ Kmat[1:(ny - 1), :] .+ 1 ./ Kmat[2:ny, :])

    TLmat = fill(NaN, ny, nx)
    # left boundary (Dirichlet)
    TLmat[:, 1] .= (2*dy/dx) .* Kmat[:, 1]
    TLmat[:, 2:nx] .= (2*dy/dx) .* K_har_x

    TRmat = fill(NaN, ny, nx)
    # right boundary (Dirichlet)
    TRmat[:, nx] .= (2*dy/dx) .* Kmat[:, nx]
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

# advection-dispersion with 1st-order upwinding (Euler)
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

    # determine dt
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

### ONLY advection, 1st-order upwind, FE
# C0: initial concentration field; c_in: inflow concentration on the left boundary.
# tedges: times [s] at which the run is split into bins; the outlet mass leaving
#   this is for the suggestion from P4-3 to use the same-time bin when comparing.
#   During each bin is accumulated so the BTC can be binned like a histogram.
# tsnap: the subset of tedges at which to store the concentration field.
function solve_advection(flow, dx, dy, ϕ, C0, c_in, tedges, tsnap)
    (; UL, UR, UU, UD) = flow
    ny, nx = size(UL)

    Ux = 0.5 .* (UL .+ UR)
    Uy = 0.5 .* (UU .+ UD)
    V_mag = sqrt.(Ux .^ 2 .+ Uy .^ 2) ./ ϕ
    dt_max = 0.1 * min(dx, dy) / maximum(V_mag)

    Cmat = copy(C0)
    Cnew = similar(Cmat)
    t = 0.0
    Cmat_ls = Matrix{Float64}[]
    Mout = zeros(length(tedges) - 1) # mass leaving the outlet in each bin
    for (n, ts) in enumerate(tedges)
        while t < ts
            dt = min(dt_max, ts - t)
            for j in 1:nx, i in 1:ny
                c = Cmat[i, j]
                # neighbours: c_in at the inlet, zero gradient at the outlet
                cL = j == 1 ? c_in : Cmat[i, j - 1]
                cR = j == nx ? c : Cmat[i, j + 1]
                cU = i == 1 ? 0.0 : Cmat[i - 1, j]
                cD = i == ny ? 0.0 : Cmat[i + 1, j]

                # upwind advective fluxes (positive = +x or +y)
                FL = UL[i, j] * (UL[i, j] >= 0 ? cL : c)
                FR = UR[i, j] * (UR[i, j] >= 0 ? c : cR)
                FU = UU[i, j] * (UU[i, j] >= 0 ? cU : c)
                FD = UD[i, j] * (UD[i, j] >= 0 ? c : cD)

                netF = -FL * dy - FU * dx + FD * dx + FR * dy
                Cnew[i, j] = c - dt / (ϕ * dx * dy) * netF

                # solute mass through the outlet face during this step
                j == nx && (Mout[n - 1] += FR * dy * dt)
            end
            Cmat, Cnew = Cnew, Cmat
            t += dt
        end
        any(s -> isapprox(s, ts), tsnap) && push!(Cmat_ls, copy(Cmat))
    end
    return (; Cmat_ls, Mout, V_mag, dt_max)
end

# Pollock method:
# v varies linearly in each direction: v = v0 + A (x - x0),
# => x(s) - x0 = v0 (exp(A s) - 1) / A 
# => time = log(v_exit / v0) / A
# xp0, yp0: start positions [m] (y = 0 on the top boundary)
# tsnap: [s] at which positions are recorded (NaN once a particle left).
# Returns the snapshots and each particle outlet arrival time
function track_particles(flow, dx, dy, ϕ, xp0, yp0, tsnap)
    (; UL, UR, UU, UD) = flow
    ny, nx = size(UL)
    np = length(xp0)

    xsnap = fill(NaN, np, length(tsnap))
    ysnap = fill(NaN, np, length(tsnap))
    tarrive = fill(Inf, np)

    # displacement after time s, starting with velocity v and gradient A
    move(v, A, s) = abs(A * s) < 1e-12 ? v * s : v * expm1(A * s) / A

    # time to reach the downstream face at signed distance d
    function exit_time(v, A, d)
        v == 0 && return Inf
        r = A * d / v
        r <= -1 && return Inf # velocity changes sign before the face
        return abs(r) < 1e-12 ? d / v : log1p(r) / A
    end

    for p in 1:np
        x, y, t = xp0[p], yp0[p], 0.0
        ix = clamp(floor(Int, x / dx) + 1, 1, nx)
        iy = clamp(floor(Int, y / dy) + 1, 1, ny)
        k = 1 # next snapshot
        while true
            xL = (ix - 1) * dx
            yU = (iy - 1) * dy
            vL, vR = UL[iy, ix] / ϕ, UR[iy, ix] / ϕ
            vU, vD = UU[iy, ix] / ϕ, UD[iy, ix] / ϕ

            Ax = (vR - vL) / dx
            Ay = (vD - vU) / dy
            vx = vL + Ax * (x - xL)
            vy = vU + Ay * (y - yU)

            # time to exit the cell through the left or right face
            τx = exit_time(vx, Ax, vx > 0 ? xL + dx - x : xL - x)
            τy = exit_time(vy, Ay, vy > 0 ? yU + dy - y : yU - y)
            τ = min(τx, τy)

            # record snapshot positions passed during this cell visit
            while k <= length(tsnap) && tsnap[k] <= t + τ
                s = tsnap[k] - t
                xsnap[p, k] = x + move(vx, Ax, s)
                ysnap[p, k] = y + move(vy, Ay, s)
                k += 1
            end
            isinf(τ) && break # stagnant: never reaches a face

            # move to the exit face and step into the neighbouring cell
            if τx <= τy
                x = vx > 0 ? xL + dx : xL
                y = clamp(y + move(vy, Ay, τ), yU, yU + dy)
                ix += vx > 0 ? 1 : -1
            else
                y = vy > 0 ? yU + dy : yU
                x = clamp(x + move(vx, Ax, τ), xL, xL + dx)
                iy += vy > 0 ? 1 : -1
            end
            t += τ

            if ix > nx # left through the outlet
                tarrive[p] = t
                break
            end
            (ix < 1 || iy < 1 || iy > ny) && break # cannot happen with these BCs
        end
    end
    return (; xsnap, ysnap, tarrive)
end
