% GPIS learning
% 2D example, Sphere mean and normals

close all; clear all; clc;

X = [-0.5,-0.5;
     -0.5,0.5;
     0.5,0.5]';
 
f = [   0,-cos(45/180*pi),-sin(45/180*pi),...
        0, -cos(45/180*pi), sin(45/180*pi),...
        0, cos(45/180*pi), sin(45/180*pi)]';

data = reshape(f, [3,length(X)])';

m = length(X); 

[Xg,Yg] = meshgrid(-1.4:0.2:1.4,-1.4:0.2:1.4);
[d1,d2] = size(Xg);
Xs = [reshape(Xg,d1*d2,1),reshape(Yg,d1*d2,1)]';
n = length(Xs);

sigma = 1;
gamma = 1;

display('Computing covariance matrix K');
K = ComputeFullKder(sigma, gamma, X, 0.2, 0);

display('Computing covariance matrix Ks');
Ks = ComputeKderX1X2(sigma, gamma, Xs, X);

display('Computing covariance matrix Kss');
Kss = ComputeFullKder(sigma, gamma, Xs, 0.2, 0.0);

figure();
imagesc(K);
colorbar;

figure();
imagesc(Ks);
colorbar;

figure();
imagesc(Kss);
colorbar;

Ks = Ks';


display('Computing means');

R = sqrt(0.5^2 + 0.5^2);
cen = [0.0, 0.0]';
mean = @(x) 1/2/R*((x-cen)'*(x-cen) - R^2);
meandx = @(x) 1/2/R*(2*(x(1)-cen(1)));
meandy = @(x) 1/2/R*(2*(x(2)-cen(2)));

for i = 1:m
    mu((i-1)*3 +1) = mean(X(:,i));
    mu((i-1)*3 +2) = meandx(X(:,i));
    mu((i-1)*3 +3) = meandy(X(:,i));
end

for i = 1:n
    mus((i-1)*3 +1) = mean(Xs(:,i));
    mus((i-1)*3 +2) = meandx(Xs(:,i));
    mus((i-1)*3 +3) = meandy(Xs(:,i));
end

mu = mu';
mus = mus';


display('Computing Regression');
fs = mus + Ks'*inv(K)*(f - mu);
sig = Kss' - Ks'*inv(K)*Ks;

figure();
imagesc(sig);
colorbar;

d = diag(sig);
devFsP = fs + d;
devFsP = reshape(devFsP(1:3:end),d1,d2);

devFsN = fs - d;
devFsN = reshape(devFsN(1:3:end),d1,d2);

Fs = reshape(fs(1:3:end),d1,d2);

figure();
hold on;
plot(X(1,:), X(2,:), 'r.', 'MarkerSize',40);
% surf(Xg,Yg,Fs);
contour(Xg,Yg,Fs,[0 0], 'LineWidth',2,'color', 'r');
quiver(X(1,:)', X(2,:)', data(:,2), data(:,3));


figure();
hold on;
plot(X(1,:), X(2,:), 'r.', 'MarkerSize',40);
surf(Xg,Yg,Fs);
contour(Xg,Yg,Fs,[0 0], 'LineWidth',2,'color', 'r');

figure();
hold on;
plot(X(1,:), X(2,:), 'r.', 'MarkerSize',40);
surf(Xg,Yg,Fs);
contour(Xg,Yg,Fs,[0 0], 'LineWidth',2,'color', 'r');
surf(Xg,Yg,devFsN);
surf(Xg,Yg,devFsP);

% Compute the inside probability.
cdfBell = @(x) 0.5.*(1 + sign(x).*sqrt(1 - exp(-2/pi.*x.*x)));

D = reshape(d(1:3:end), d1,d1);
prob = cdfBell((0-Fs)./D);
figure();
hold on;
surf(Xg,Yg,prob);

