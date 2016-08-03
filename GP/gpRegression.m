close all; clear all; clc;

l = 0.3;
kernel = @(x,y) 0.05*exp(-(x-y)*(x-y)'/(l*l));

%% Let be X the finite set of variables that we want to sample.
X = 0:0.01:1.5
n = length(X);

K = [];
for i=1:n
    for j =1:n
        K(i,j) = kernel(X(i),X(j));
    end
end

fs = [0, -0.3, 0];
Xs = [0.4, 0.6, 0.9];

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

figure(2);
hold on;
plot(X, nus,'b','LineWidth',2);
plot(Xs,fs, 'r*');

d = diag(sig);
plot(X, nus+sqrt(d),'k');
plot(X, nus-sqrt(d),'k');

ix1 = 10;
x1 = X(ix1);
nu1 = nus(ix1);
sig1 = sig(ix1,ix1);
f1 = -0.6:0.01:0.6;
g1 = normpdf(f1, nu1, sqrt(sig1))*0.1;
plot(x1 + g1, f1,'r-')
plot(x1*ones(1, length(f1)), f1, 'k--');

ix2 = 80;
x2 = X(ix2);
nu2 = nus(ix2);
sig2 = sig(ix2,ix2);
f2 = -0.6:0.01:0.6;
g2 = normpdf(f2, nu2, sqrt(sig2))*0.1;
plot(x2 + g2, f2,'r-')
plot(x2*ones(1, length(f2)), f2, 'k--');

% Horizontal line
plot(X, zeros(1, length(X)), 'k-');