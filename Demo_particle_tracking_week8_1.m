clc; clf; clear all;
%% physical properties
rho = 1000; % water density [kg/m^3]
mu = 1e-3; % dynamic viscosity [Pa s]
g = 9.81; % gravitational acceleration [m/s^2]
k = 1e-10; % intrinsic permeability [m^2]

%% model domain
Lx = 100; % [m]
Ly = 50; % [m]

nx = 100; 
ny = 50;
nc = nx * ny;

dx = Lx/nx;
dy = Ly/ny;


%% boundary  conditions
uin = 0.1 / 100 / 60; % [m/s] == 0.1 cm/min   <-- Neumann on left boundary
hR = 9; % [m] <-- Dirichelet on right boundary

%% Transmissibility matrix
clc
% Kmat = rho*g*k/mu * ones(ny,nx);
kmat = readmatrix('permeability_field.csv');


kmean = mean(log10(kmat(:)));
kstd = std(log10(kmat(:)));

Y = (log10(kmat) - kmean)/kstd;

kmat = kmean + 1*Y;

Kmat = rho*g/mu * 10.^kmat;

K_har_x = (1./Kmat(:,1:nx-1) + 1./Kmat(:,2:nx)).^(-1);
K_har_y = (1./Kmat(1:ny-1,:) + 1./Kmat(2:ny,:)).^(-1);

TLmat = nan(ny,nx);
TLmat(:,1) = 0;
TLmat(:,2:nx) = 2/dx*dy * K_har_x;

TRmat = nan(ny,nx);
TRmat(:,nx) = 2/dx*dy*Kmat(:,nx);
TRmat(:,1:nx-1) = 2/dx*dy * K_har_x;

TUmat = nan(ny,nx);
TUmat(1,:) = 0;
TUmat(2:ny,:) = 2/dy*dx * K_har_y;

TDmat = nan(ny,nx);
TDmat(ny,:) = 0;
TDmat(1:ny-1,:) = 2/dy*dx * K_har_y;

TL = TLmat(:);
TR = TRmat(:);
TU = TUmat(:);
TD = TDmat(:);

TC = TL + TR + TU + TD;

Tmat = diag(TC,0) + diag(-TL(ny+1:nc),-ny) + diag(-TU(2:nc),-1)...
    +diag(-TD(1:nc-1),1) + diag(-TR(1:nc-ny),ny);

%% RHS vector b
bmat = zeros(ny,nx);
bmat(:,1) = uin*dy;
bmat(:,nx) = TRmat(:,nx)*hR;

b = bmat(:);

%% head vector h
h = Tmat\b;

hmat = reshape(h,ny,nx);


subplot(211)
imagesc(log10(Kmat))
colorbar

subplot(212)
imagesc(hmat)
colorbar

%% Darcy velocity

dh_x = hmat(:,1:nx-1)-hmat(:,2:nx);
dh_y = hmat(1:ny-1,:)-hmat(2:ny,:);

UL = nan(ny,nx);
UL(:,1) = uin;
UL(:,2:nx) = TLmat(:,2:nx)/dy.*dh_x;

UR = nan(ny,nx);
UR(:,nx) = TRmat(:,nx)/dy .* (hmat(:,nx) - hR);
UR(:,1:nx-1) = TRmat(:,1:nx-1)/dy .* dh_x;

UU = nan(ny,nx);
UU(1,:) = 0;
UU(2:ny,:) = TUmat(2:ny,:)/dx .* dh_y;

UD = nan(ny,nx);
UD(ny,:) = 0;
UD(1:ny-1,:) = TDmat(1:ny-1,:)/dx .* dh_y;

Ux = 0.5 * (UL + UR);
Uy = 0.5 * (UU + UD);

%% visualization


subplot(211)
imagesc(log10(Kmat)); colorbar;

axis equal tight

subplot(212)
imagesc(hmat); colorbar;


hold on;

% Plot only every few velocity vectors
xCells = dx/2:dx:Lx;
yCells = dy/2:dy:Ly;

[xx,yy] = meshgrid(xCells,yCells);

skip = 5;

quiver( ...
    xx(1:skip:end,1:skip:end), ...
    yy(1:skip:end,1:skip:end), ...
    Ux(1:skip:end,1:skip:end), ...
    Uy(1:skip:end,1:skip:end), ...
    'k');

xlabel('x [m]');
ylabel('y [m]');
axis equal tight;

%% particle tracking
clc
% intial setting
phi = 0.3; % porosity

np = 20; % number of particles

step = 1; % iteration step
maxStep = 200;

xTrack = nan(np,maxStep);
yTrack = nan(np,maxStep);
tTrack = nan(np,maxStep);

xTrack(:,1) = zeros(np,1)+1e-10*dx;
yTrack(:,1) = linspace(dy/10, Ly-dy/10,np);
tTrack(:,1) = 0;

active = true(np,1); % whether still inside the model domain
while any(active)
%% next location
    xp = xTrack(:,step);
    yp = yTrack(:,step);
    tp = tTrack(:,step);
    
    ip = find(active);
    
    xp = xp(ip);
    yp = yp(ip);
    tp = tp(ip);
    
    % identify index of the cells containing particles
    ix = floor(xp/dx)+1;
    iy = floor(yp/dy)+1;
    
    % identify cell face pore velocities
    ind = ny*(ix-1)+iy;
    vL = UL(ind)/phi;
    vU = UU(ind)/phi; % pore velocity on Up face
    vD = UD(ind)/phi;
    vR = UR(ind)/phi;
    
    % identify cell face coordinates
    xL = (ix-1)*dx;
    xR = ix*dx;
    yU = (iy-1)*dy;
    yD = iy*dy;
    
    % velociy gradiets
    Ax = (vR - vL) / dx;
    Ay = (vD - vU) / dy;
    
    % velocity on particle locations
    vxp = vL + Ax.*(xp - xL);
    vyp = vU + Ay.*(yp - yU);
    
    % exit velocity
    vxExit = zeros(size(xp));
    vyExit = zeros(size(yp));
    
    
    ind = vxp >= 0;% in x direction
    vxExit(ind) = vR(ind);
    vxExit(~ind) = vL(~ind);
    
    ind = vyp >= 0;
    vyExit(ind) = vD(ind);
    vyExit(~ind) = vU(~ind);
    
    % exit time
    tauX = 1./Ax .* log(vxExit ./ vxp);
    tauY = 1./Ay .* log(vyExit ./ vyp);
    tau = min(tauX, tauY);
    
    % move particle
    xp = xL + (vxp .* exp(Ax.*tau) - vL) ./ Ax;
    yp = yU + (vyp .* exp(Ay.*tau) - vU) ./ Ay;
    tp = tp + tau;
    
    % infinitesimally move particle into next cell (nudging)
    ind = vxp >= 0;% in x direction
    xp(ind) = xp(ind) + 1e-10*dx;
    xp(~ind) = xp(~ind) - 1e-10*dx;
    
    
    ind = vyp >= 0;% in x direction
    yp(ind) = yp(ind) + 1e-10*dy;
    yp(~ind) = yp(~ind) - 1e-10*dy;
    
    % save
    step = step +1;
    xTrack(ip,step) = xp;
    yTrack(ip,step) = yp;
    tTrack(ip,step) = tp;
    
    % check if particles reached the outlet
    ind = xTrack(ip,step) > Lx - 1e-10*dx;
    active(ip(ind)) = false;

end
% visualization
clf

subplot(211)
imagesc(log10(Kmat)); colorbar;

hold on
for p = 1:np
    plot(xTrack(p,:),yTrack(p,:),'k-','linewidth',1.2)
end

axis equal tight

subplot(212)
imagesc(hmat); colorbar;


hold on;

% Plot only every few velocity vectors
xCells = dx/2:dx:Lx;
yCells = dy/2:dy:Ly;

[xx,yy] = meshgrid(xCells,yCells);

skip = 5;

quiver( ...
    xx(1:skip:end,1:skip:end), ...
    yy(1:skip:end,1:skip:end), ...
    Ux(1:skip:end,1:skip:end), ...
    Uy(1:skip:end,1:skip:end), ...
    'k');

xlabel('x [m]');
ylabel('y [m]');
axis equal tight;
drawnow
