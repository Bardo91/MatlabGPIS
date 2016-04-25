%% Learning GP
% GP sampling
close all; clear all; clc;

l = 3;
kernel = @(x,y) exp(-2*sin(5*pi*(x-y))/l/l);

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

figure(1);
imagesc(K);
colorbar;


%% Sample gaussians
figure(2);
hold on;
for i=1:5
    u = randn(n,1);
    % L  = chol(K); not very stable numerically 
    [A S D] = svd(K);
    L = A*sqrt(S);
    F = L*u;
    plot(F);
end
