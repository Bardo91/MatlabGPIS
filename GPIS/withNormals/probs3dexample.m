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

X = appleLoc(1:50:end,:)';
norms = appleNorm(1:50:end,:)';

% X = X(:,X(1,:) < 0);
% norms = norms(:,X(1,:) < 0);

absmax = max(abs(X'))';
m = length(X); 
for i = 1:m
    X(:,i) = X(:,i)./absmax;
end

f = zeros(m,1);
f = [f,norms'];
[a b] = size(f);
f = reshape(f', a*b, 1);

step = 0.25;
lim=2;
%  spherical mean
% sigma = 0.0098;
% gamma = 400;    
% noiseVal = 5e-6;
% noiseGrad = 0.02;
% R = 1;
% cen = [0.0, 0.0, 0.0]';
% mean = @(x)1/2/R*((x-cen)'*(x-cen) - R^2);
% meandx = @(x) 1/R*((x(1)-cen(1)));
% meandy = @(x) 1/R*((x(2)-cen(2)));
% meandz = @(x) 1/R*((x(3)-cen(3)));

% No prior
sigma = 0.75;%0.5395;
gamma = 10;%10.43;    
noiseVals = 1e-3;
noiseGrad = 0.03;
meanLevel = 0.1;

mean = @(x)meanLevel;
meandx = @(x)0;
meandy = @(x)0;
meandz = @(x)0;

[Xg,Yg, Zg] = meshgrid(-lim:step:lim,-lim:step:lim,-lim:step:lim);
[d1,d2] = size(Xg);
Xs = [reshape(Xg,d1*d2,1),reshape(Yg,d1*d2,1),reshape(Zg,d1*d2,1)]';
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
sig = Kss - Ks*kinv*Ks';


Fs = reshape(fs(1:4:end),d1,d1,d1);

d = diag(sig);
devFsP = fs + d;
devFsP = reshape(devFsP(1:4:end),d1,d1,d1);

devFsN = fs - d;
devFsN = reshape(devFsN(1:4:end),d1,d1,d1);

figure();
hold on;
plot3(X(1,:), X(2,:), X(3,:), 'r.', 'MarkerSize',40);
quiver3(X(1,:), X(2,:), X(3,:), norms(1,:), norms(2,:), norms(3,:));

p = patch(isosurface(Xg, Yg, Zg, Fs, 0));
p.FaceColor = 'green';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3)
camlight; lighting phong;



display('Computing the inside probability');
% Compute the inside probability.
cdfBell = @(x) 0.5.*(1 + sign(x).*sqrt(1 - exp(-2/pi.*x.*x)));

D = reshape(d(1:4:end), d1,d1,d1);
prob = cdfBell((0-Fs)./D);
figure();
hold on;
p0 = patch(  isosurface(Xg, Yg, Zg, prob, 0), ...
                'FaceColor','cyan',...
                'FaceAlpha',0.25,...
                'EdgeColor', 'none');
p25 = patch(  isosurface(Xg, Yg, Zg, prob, 0.25), ...
                'FaceColor','blue',...
                'FaceAlpha',0.25,...
                'EdgeColor', 'none');
p50 = patch(  isosurface(Xg, Yg, Zg, prob, 0.5), ...
                'FaceColor','green',...
                'FaceAlpha',0.25,...
                'EdgeColor', 'none');
p75 = patch(  isosurface(Xg, Yg, Zg, prob, 0.75), ...
                'FaceColor','red',...
                'FaceAlpha',0.25,...
                'EdgeColor', 'none');
p1 = patch(  isosurface(Xg, Yg, Zg, prob, 1), ...
                'FaceColor','yellow',...
                'FaceAlpha',0.5,...
                'EdgeColor', 'none');
daspect([1 1 1])
plot3(X(1,:), X(2,:), X(3,:), 'r.', 'MarkerSize',40);
view(3)
camlight; lighting phong;


figure()
hold on;
plot3(X(1,:), X(2,:), X(3,:), 'r.', 'MarkerSize',40);
p = patch(isosurface(Xg, Yg, Zg, Fs, 0),'FaceAlpha',0.7);
p.FaceColor = 'green';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3)
camlight; lighting phong;
colormap('default');   % set colormap
contourf(Xg(:,:,1),Yg(:,:,1),prob(:,:,floor(d1/2)));        % draw image and scale colormap to values range
colorbar;