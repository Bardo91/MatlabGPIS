%% Learning GP
% GP regression
close all; clear all; clc;

l = 0.1;
kernel = @(x,y) exp(-(x-y)*(x-y)'/(l*l));

%% Let be X the finite set of variables that we want to sample.
X = 0:0.01:1;

% Assuming 0 mean. Compute the covariance assuming exponential kernel
n = length(X);

K = [];
for i=1:n
    for j =1:n
        K(i,j) = kernel(X(i),X(j));
    end
end

fs = [0, 1, 0, 4];
Xs = [0, 0.2, 0.6, 0.8];

% fs = [0, 1, 0, 4, 2];
% Xs = [0, 0.2, 0.6, 0.8, 0.5];


Kss = [];
for i=1:length(fs)
    for j =1:length(fs)
        Kss(i,j) = kernel(Xs(i),Xs(j));
    end
end

Ks = [];
for i=1:length(fs)
    for j =1:n
        Ks(i,j) = kernel(Xs(i),X(j));
    end
end

nus = 0 + Ks'*inv(Kss)*(fs');
sig = K - Ks'*inv(Kss)*Ks;

figure(1);
imagesc(sig);
colorbar;

figure(2);
hold on;
plot(X, nus);
plot(Xs,fs, 'r*');

d = diag(sig);
plot(X, nus+sqrt(d), 'r');
plot(X, nus-sqrt(d), 'r');

figure(3);
hold on;
for i=1:5
    u = randn(n,1);
    % L  = chol(K); not very stable numerically 
    [A S D] = svd(sig);
    L = A*sqrt(S);
    F = L*u;
    plot(nus + F);
end