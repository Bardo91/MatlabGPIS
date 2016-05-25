% GPIS learning
% 1D example, zero mean

close all; clear all; clc;
% 
% X = [1,0,0;
%     0,1,0;
%     0,0,1;
% %     -0.5774,-0.5774,-0.5774]';
%     -1,-1,-1]';
% 
% m = length(X); 
% norms = [   1 0 0 -sqrt(3);
%             0 1 0 -sqrt(3);
%             0 0 1 -sqrt(3)]';
% f = [   0,1,0,0,...
%         0,0,1,0,...
%         0,0,0,1,...
%         0,-sqrt(3),-sqrt(3),-sqrt(3)   ]';
% 
% step = 0.5;
% lim=2;
% sigma = 0.5;
% gamma = 1;
% R = 1;

AppleData;

X = appleLoc(1:20:end,:)';
norms = appleNorm(1:20:end,:)';

% norms = norms(:,X(1,:) < 0);
% X = X(:,X(1,:) < 0);

absmax = max(abs(X'))';
m = length(X); 
for i = 1:m
    X(:,i) = X(:,i)./absmax;
end

f = zeros(m,1);
f = [f,norms'];
[a b] = size(f);
f = reshape(f', a*b, 1);

step = 0.1;
lim=1.3;

% %  spherical mean
% sigma = 0.4; 
% gamma = 100; 
% noiseVals = 5e-6; 
% noiseGrad = 0.02; 
% R = 1;
% cen = [0.0, 0.0, 0.0]';
% mean = @(x)1/2/R*((x-cen)'*(x-cen) - R^2);
% meandx = @(x) 1/R*((x(1)-cen(1)));
% meandy = @(x) 1/R*((x(2)-cen(2)));
% meandz = @(x) 1/R*((x(3)-cen(3)));
% prior = struct( 'pos',[0 0 0],'type','N',...
%                 'param', [1 1 1], 'rot', [2*pi 0 0],...
%                 'Sigma', sigma, 'Gamma', gamma,...
%                 'noiseVals', noiseVals, 'noiseGrad', noiseGrad);
            
% No prior
sigma = 0.11;%0.5395;
gamma = 11.2;%10.43;    
noiseVals = 1.17e-5;
noiseGrad = 5.87e-4;
meanLevel = 0.17;

prior = struct( 'pos',[0 0 0],'type','N',...
                'param', [meanLevel 1 1], 'rot', [2*pi 0 0],...
                'Sigma', sigma, 'Gamma', gamma,...
                'noiseVals', noiseVals, 'noiseGrad', noiseGrad);

mean = @(x)meanLevel;
meandx = @(x)0;
meandy = @(x)0;
meandz = @(x)0;

[Xg,Yg,Zg] = meshgrid(-lim:step:lim,-lim:step:lim,[-0.8, 0, 0.8]);
[d1,d2,d3] = size(Xg);
Xs = [reshape(Xg,d1*d2*d3,1),reshape(Yg,d1*d2*d3,1),reshape(Zg,d1*d2*d3,1)]';
n = length(Xs);

display('Computing covariance matrix K');
K = ComputeFullKder(sigma, gamma, X, 0.0, 0.0);

display('Computing covariance matrix Ks');
Ks = ComputeKderX1X2(sigma, gamma, Xs, X);

display('Computing covariance matrix Kss');
[D,N] = size(Xs);
Kss = DiagComputeKderX1X2(sigma,gamma,Xs,Xs);
NoiseDiag =  [noiseVals * ones(1,N);
              noiseGrad * ones(D,N)];
Kss = Kss + diag(NoiseDiag(:));

display('Computing means');

for i = 1:m
    mu((i-1)*4 +1) = mean(X(:,i));
    mu((i-1)*4 +2) = meandx(X(:,i));
    mu((i-1)*4 +3) = meandy(X(:,i));
    mu((i-1)*4 +4) = meandz(X(:,i));
end

for i = 1:n
    mus((i-1)*4 +1) = mean(Xs(:,i));
    mus((i-1)*4 +2) = meandx(Xs(:,i));
    mus((i-1)*4 +3) = meandy(Xs(:,i));
    mus((i-1)*4 +4) = meandz(Xs(:,i));
end

mu = mu';
mus = mus';

display('Computing regression');
kinv = inv(K);
fs = mus + Ks*kinv*(f - mu);
sig = - Ks*kinv*Ks';

d = diag(Kss) - diag(sig);

Fs = reshape(fs(1:4:end),d1,d2,d3);

display('Computing the inside probability');
% Compute the inside probability.
cdfBell = @(x) 0.5.*(1 + sign(x).*sqrt(1 - exp(-2/pi.*x.*x)));

D = reshape(d(1:4:end), d1,d2,d3);
prob = cdfBell((0-Fs)./D);



[meanValue, meanGrad] = computePriorFunctions(prior)
%% Accurate plot
[faces, vertices] = computeSurface(X, norms, prior, meanValue, meanGrad, X(:,1:20:end), 0.05, false);

figure
hold on
axis equal

plot3(X(1,:),X(2,:),X(3,:),'r.','markersize',30);
quiver3(X(1,:),X(2,:),X(3,:), norms(1,:),norms(2,:),norms(3,:),'linewidth',2,'color','r');

patch('faces',faces,'vertices',vertices,...
    'facecolor',[0.5 0.5 0.5], ...
    'edgecolor', 'none', ...
    'facelighting','phong',...
    'FaceAlpha', 0.5);
camlight
set(gca,'view',[46.8000   18.8000]);
light('Position',[-1 -1 0])
contourf(Xg(:,:,1),Yg(:,:,1),prob(:,:,2));        % draw image and scale colormap to values range
colorbar;
view([-37.5 30])

figure
% subplot(1,4,1)
hold on
axis equal
plot3(X(1,:),X(2,:),X(3,:),'r.','markersize',30);
quiver3(X(1,:),X(2,:),X(3,:), norms(1,:),norms(2,:),norms(3,:),'linewidth',2,'color','r');
patch('faces',faces,'vertices',vertices,...
    'facecolor',[0.5 0.5 0.5], ...
    'edgecolor', 'none', ...
    'facelighting','phong',...
    'FaceAlpha', 0.5);
camlight
set(gca,'view',[46.8000   18.8000]);
light('Position',[-1 -1 0])
view([-37.5 30])

figure
% subplot(1,4,2)
contourf(Xg(:,:,1),Yg(:,:,1),prob(:,:,1));        % draw image and scale colormap to values range
colorbar;
axis equal

figure
% subplot(1,4,3)
contourf(Xg(:,:,1),Yg(:,:,1),prob(:,:,2));        % draw image and scale colormap to values range
colorbar;
axis equal

figure
% subplot(1,4,4)
contourf(Xg(:,:,1),Yg(:,:,1),prob(:,:,3));        % draw image and scale colormap to values range
colorbar;
axis equal