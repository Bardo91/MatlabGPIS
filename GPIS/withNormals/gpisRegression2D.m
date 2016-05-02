% GPIS learning
% 2D example, Sphere mean and normals

close all; clear all; clc;

X = [0,0;
     1,0;
     1,1;
     0,1]';

m = length(X); 

[Xg,Yg] = meshgrid(-2:0.25:2,-2:0.25:2);
[d1,d2] = size(Xg);
Xs = [reshape(Xg,d1*d2,1),reshape(Yg,d1*d2,1)]';
n = length(Xs);

sigma = 1;
gamma = 1;
K = ComputeFullKder(sigma, gamma, X, 0, 0);
Ks = ComputeKderX1X2(sigma, gamma, Xs, X);
Kss = ComputeKderX1X2(sigma, gamma, Xs, Xs);

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

f = [0,-1,-1,  0, 1, -1, 0, 1, 1, 0, -1, 1]';
fs = mus + Ks'*inv(K)*(f - mu);
sig = Kss' - Ks'*inv(K)*Ks;

figure();
imagesc(sig);
colorbar;

figure();
hold on;
plot(X(1,:), X(2,:), 'r*');

d = diag(sig);
devFsP = fs + d;
devFsP = reshape(devFsP(1:3:end),d1,d2);

devFsN = fs - d;
devFsN = reshape(devFsN(1:3:end),d1,d2);

Fs = reshape(fs(1:3:end),d1,d2);
surf(Xg,Yg,Fs);
surf(Xg,Yg,devFsN);
surf(Xg,Yg,devFsP);
axis([-2 0.5 -2 2 -2 2]);

% Compute the inside probability.
cdfBell = @(x) 0.5.*(1 + sign(x).*sqrt(1 - exp(-2/pi.*x.*x)));

D = reshape(d(1:3:end), d1,d1);
prob = cdfBell((0-Fs)./D);
figure();
hold on;
surf(Xg,Yg,prob);

