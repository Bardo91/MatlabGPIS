% GPIS learning
% 1D example, zero mean

% close all; clear all; clc;

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

f = [0,-1,-1,  0, 1, -1, 0, 1, 1, 0, -1, 1]';
fs = 0 + Ks'*inv(K)*(f);
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
