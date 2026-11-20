function [hmat] = forwardModel(kmat)
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
% clc
% % Kmat = rho*g*k/mu * ones(ny,nx);
% kmat = readmatrix('permeability_field.csv');
% 
% 
% kmean = mean(log10(kmat(:)));
% kstd = std(log10(kmat(:)));
% 
% Y = (log10(kmat) - kmean)/kstd;
% 
% kmat = kmean + 0.1*Y;
% 
% Kmat = rho*g/mu * 10.^kmat; % hydraulic conductivity

Kmat = rho*g*kmat/mu;

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

Tmat = sparse(Tmat);
%% RHS vector b
bmat = zeros(ny,nx);
bmat(:,1) = uin*dy;
bmat(:,nx) = TRmat(:,nx)*hR;

b = bmat(:);

%% head vector h
h = Tmat\b;

hmat = reshape(h,ny,nx);


% subplot(211)
% imagesc(log10(Kmat))
% colorbar
% 
% subplot(212)
% imagesc(hmat)
% colorbar

%% Darcy velocity

% dh_x = hmat(:,1:nx-1)-hmat(:,2:nx);
% dh_y = hmat(1:ny-1,:)-hmat(2:ny,:);
% 
% UL = nan(ny,nx);
% UL(:,1) = uin;
% UL(:,2:nx) = TLmat(:,2:nx)/dy.*dh_x;
% 
% UR = nan(ny,nx);
% UR(:,nx) = TRmat(:,nx)/dy .* (hmat(:,nx) - hR);
% UR(:,1:nx-1) = TRmat(:,1:nx-1)/dy .* dh_x;
% 
% UU = nan(ny,nx);
% UU(1,:) = 0;
% UU(2:ny,:) = TUmat(2:ny,:)/dx .* dh_y;
% 
% UD = nan(ny,nx);
% UD(ny,:) = 0;
% UD(1:ny-1,:) = TDmat(1:ny-1,:)/dx .* dh_y;
% 
% Ux = 0.5 * (UL + UR);
% Uy = 0.5 * (UU + UD);

%% visualization


% subplot(211)
% imagesc(log10(Kmat)); colorbar;
% 
% axis equal tight
% 
% subplot(212)
% imagesc(hmat); colorbar;
% 
% 
% hold on;
% 
% % Plot only every few velocity vectors
% xCells = dx/2:dx:Lx;
% yCells = dy/2:dy:Ly;
% 
% [xx,yy] = meshgrid(xCells,yCells);
% 
% skip = 5;
% 
% quiver( ...
%     xx(1:skip:end,1:skip:end), ...
%     yy(1:skip:end,1:skip:end), ...
%     Ux(1:skip:end,1:skip:end), ...
%     Uy(1:skip:end,1:skip:end), ...
%     'k');
% 
% xlabel('x [m]');
% ylabel('y [m]');
% axis equal tight;

end