clc; clf; clear all;
%% physical properties
rho = 1000; % water density [kg/m^3]
mu = 1e-3; % dynamic viscosity [Pa s]
g = 9.81; % gravitational acceleration [m/s^2]
k = 1e-10; % intrinsic permeability [m^2]


phi = 0.3; % porosity [-]
D_diff = 1.e-9; % molecular diffusion coeff. [m^2/s]
alpha = 1; % isotropiic dispersivity [m]
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

kmat = kmean + 0.5*Y;

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

%% pore velocity
V_mag = sqrt(Ux.^2 + Uy.^2) / phi;% field pore-velocity magnitude
V_max = max(V_mag(:));

%% transport (advection + dispersion)
clc

t = 0; % time [s]
tEndDay = 200;

Cmat = zeros(ny,nx);

C_in = 1; % left boundary concentration (Dirichlet)

doPlot = false; % should I plot?
tPlot = [1:tEndDay]*24*60*60;% when to plot
itPlot = 1; % counting plots

Mout = nan(length(tPlot),1);

while t < tEndDay*24*60*60
    if t > 5*24*60*60
        C_in = 0;
    end

    
    %% neighboring cell centers' concentration
    CL = nan(ny,nx);
    CR = nan(ny,nx);
    CU = nan(ny,nx);
    CD = nan(ny,nx);
    
    CL(:,1) = ones(ny,1)*C_in;
    CL(:,2:nx) = Cmat(:,1:nx-1);
    
    CR(:,nx) = Cmat(:,nx); % zero-gradient concentration condition [Neumann]
    CR(:,1:nx-1) = Cmat(:,2:nx);
    
    CU(1,:) = 0;
    CU(2:ny,:) = Cmat(1:ny-1,:);
    
    CD(ny,:) = 0;
    CD(1:ny-1,:) = Cmat(2:ny,:);
    
    %% advective flux estimation (on cell faces)
    Fadv_L = nan(ny,nx);
    Fadv_U = nan(ny,nx);
    Fadv_D = nan(ny,nx);
    Fadv_R = nan(ny,nx);
    
    % upwinding
    % left 
    ind = UL >= 0;
    Fadv_L(ind) = UL(ind).*CL(ind);
    Fadv_L(~ind) = UL(~ind).*Cmat(~ind);
    
    % right
    ind = UR >= 0;
    Fadv_R(ind) = UR(ind).*Cmat(ind);
    Fadv_R(~ind) = UR(~ind).*CR(~ind);
    
    % up
    ind = UU >= 0;
    Fadv_U(ind) = UU(ind).*CU(ind);
    Fadv_U(~ind) = UU(~ind).*Cmat(~ind);
    
    % down
    ind = UD >= 0;
    Fadv_D(ind) = UD(ind).*Cmat(ind);
    Fadv_D(~ind) = UD(~ind).*CD(~ind);
    
    % net flux
    netFadv = -Fadv_L*dy -Fadv_U*dx +Fadv_D*dx +Fadv_R*dy;

    %% dispersive flux estimation (on cell faces)
    % clc
    % 
    % % dispersion coeff on cell centers (D: D_diff + alpha * V_mag)
    % Dmat = D_diff + alpha * V_mag;
    % 
    % % dispersion coeff on internal faces (exluding domain boundary faces)
    % Dx_har = 2./(1./Dmat(:,1:nx-1) +  1./Dmat(:,2:nx));
    % Dy_har = 2./(1./Dmat(1:ny-1,:) +  1./Dmat(2:ny,:));
    % 
    % % dispsersion coeff on cell faces
    % DL = nan(ny,nx);
    % DR = nan(ny,nx);
    % DU = nan(ny,nx);
    % DD = nan(ny,nx);
    % 
    % DL(:,1) = 2 * Dmat(:,1);
    % DL(:,2:nx) = Dx_har;
    % 
    % DR(:,nx) = Dmat(:,nx);
    % DR(:,1:nx-1) = Dx_har;
    % 
    % DU(1,:) = 0;
    % DU(2:ny,:) = Dy_har;
    % 
    % DD(ny,:) = 0;
    % DD(1:ny-1,:) = Dy_har;
    % 
    % % flux
    % Fdis_L = phi * DL .* (CL-Cmat)/dx;
    % Fdis_R = phi * DR .* (Cmat-CR)/dx;
    % Fdis_U = phi * DU .* (CU-Cmat)/dy;
    % Fdis_D = phi * DD .* (Cmat-CD)/dy;
    % 
    % 
    % % net flux
    % netFdis = -Fdis_L*dy -Fdis_U*dx +Fdis_D*dx +Fdis_R*dy;
    % 


    
    %% determine dt
    
    dt_adv = 0.1* dx/V_max;
    % dt_dis = 0.1* dx^2/max(Dmat(:));
    
    dt = dt_adv;%min(dt_adv,dt_dis);

    if t+dt > tPlot(itPlot)
        dt = tPlot(itPlot)-t;
        doPlot = true;
        itPlot = itPlot +1;

    end


    %% advance concentration field in time
    t = t+dt;

    
    Cmat = Cmat - dt/(phi*dx*dy) * (netFadv);%+netFdis);
    

    if doPlot == true
    
        Mout(itPlot-1) = sum(UR(:,end).*Cmat(:,end)*dy);
        %% visualization
    
        % clf
        % subplot(311)
        % imagesc(log10(Kmat)); %colorbar;
        % 
        % axis equal tight
        % 
        % subplot(312)
        % imagesc(hmat);% colorbar;
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
        % axis equal tight;
        % 
        % subplot(313)
        % imagesc(Cmat); %colorbar;
        % caxis([0 1])
        % 
        % title(sprintf('times %.2f days',t/24/60/60))
        % 
        % xlabel('x [m]');
        % ylabel('y [m]');
        % axis equal tight;
        % 
        % drawnow
        fprintf('times %.2f days\n',t/24/60/60)
        doPlot = false;
    end
end

%%
clf
Mexit = trapz(tPlot,Mout);
f_FVM = Mout/Mexit;
trapz(tPlot,f_FVM)

% plot(tPlot,f_FVM);
% xlabel('time [s]')
% ylabel('pdf [1/s]')

plot(tPlot/24/60/60,f_FVM*24*60*60);
xlabel('time [d]')
ylabel('pdf [1/d]')