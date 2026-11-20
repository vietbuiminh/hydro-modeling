clc; clear all; clf
%%
kmat_true = readmatrix('permeability_field.csv');
kmat_ensemble = readmatrix('permeability_ensemble.csv');

ny = 50;
nx = 100;

hSig = 0.02;


%% monitoring system setting
xInd = round(linspace(10,nx-9,10));
yInd = round(linspace(5,ny-4,5));

[xx,yy] = meshgrid(xInd,yInd);

obsInd = (xx(:)-1)*ny+yy(:);

Nd = length(obsInd);

%% observation data
hmat_true = forwardModel(kmat_true);
obs = hmat_true(obsInd) + hSig*randn(Nd,1);

%% predictions
yPrior = zeros(Nd,Ne);
for i=1:Ne
    kmat_member = reshape(kmat_ensemble(:,i),ny,nx);
    hmat_predict = forwardModel(kmat_member);
    yPrior(:,i) = hmat_predict(obsInd);
    fprintf('%d\n',i)
end
%% ensemble smoothiing
xPrior = log10(kmat_ensemble);
xPriorMean = mean(xPrior,2);
yPriorMean = mean(yPrior,2);

Ax = xPrior - xPriorMean;
Ay = yPrior - yPriorMean;

Cxy = (Ax*Ay') / (Ne-1);
Cyy = (Ay*Ay') / (Ne-1);
Cd = eye(Nd)*hSig^2;

K_gain = Cxy / (Cd+Cyy);

xPost = xPrior + K_gain * (obs+rand(Nd,Ne)- yPrior);


%%
Ne = size(kmat_ensemble,2);

% for i=1:Ne
    clf
    subplot(311)
    imagesc(log10(kmat_true));
    axis equal tight
    cl = caxis;
    hold on
    plot(xx(:),yy(:),'wo')
    
    subplot(312)
    kmat_mean_prior = reshape(mean(log10(kmat_ensemble),2),ny,nx);
    imagesc((kmat_mean_prior))
    axis equal tight
    caxis(cl)

    subplot(313)
    kmat_mean_post = reshape(mean(xPost,2),ny,nx);
    imagesc(kmat_mean_post)
    axis equal tight
    caxis(cl)

    

    pause(1)
% end

%% 