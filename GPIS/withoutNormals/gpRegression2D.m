% GPIS learning
% 1D example, zero mean

close all; clear all; clc;

X = [0,0;
     1,0;
     1,1;
     0,1;
     0.5,1]';
f =[0,0,0,0,0]';

m = length(X); 

[Xg,Yg] = meshgrid(-2:0.1:2,-2:0.1:2);
[d1,d2] = size(Xg);
Xs = [reshape(Xg,d1*d2,1),reshape(Yg,d1*d2,1)]';
n = length(Xs);

sigma = 1;
gamma =10;

kernel = @(x,y)(sigma^2 * exp(-1/2 * gamma *(x - y)'*(x - y)));

R = 1;
cen = [0.5, 0.5]';
mean = @(x) 1/2/R*((x-cen)'*(x-cen) - R*2);


K = [];
for i = 1:m
    for j = 1:m
       K(i,j) = kernel(X(:,i), X(:,j));
    end
end


Ks = [];
for i = 1:m
    for j = 1:n
       Ks(i,j) = kernel(X(:,i), Xs(:,j)); 
    end
end


Kss = [];
for i = 1:n
    for j = 1:n
       Kss(i,j) = kernel(Xs(:,i), Xs(:,j));
    end
end


figure();
imagesc(K);
colorbar;

figure();
imagesc(Ks);
colorbar;

figure();
imagesc(Kss);
colorbar;


for i = 1:m
    mu(i) = mean(X(:,i));
end
fss=[];
for i = 1:n
    mus(i) = mean(Xs(:,i));
end
mu = mu';
mus = mus';

fs = mus + Ks'*inv(K)*(f - mu);
sig = Kss' - Ks'*inv(K)*Ks;

figure();
imagesc(sig);
colorbar;


d = diag(sig);
devFsP = fs + d;
devFsP = reshape(devFsP,d1,d2);

devFsN = fs - d;
devFsN = reshape(devFsN,d1,d2);

Fs = reshape(fs,d1,d2);
figure();
hold on;
plot(X(1,:), X(2,:), 'r*');
surf(Xg,Yg,Fs);
surf(Xg,Yg,devFsN);
surf(Xg,Yg,devFsP);
shading interp;
