% GPIS learning
% 2D example, Sphere mean wo normals

close all; clear all; clc;

X = [-0.5,-0.5;
     -0.5,0.5;
     0.5,0.5;
%      0.5,-.5]';
     0.3,1;
     00.5,-1]';
%      0.0,0.707;
%      0.707,-0.0;
%      -0.707,0.0;
%      0.0, -0.707]';
 
m = length(X); 

f = zeros(m,1);

[Xg,Yg] = meshgrid(-1.5:0.1:1.5,-1.5:0.1:1.5);
[d1,d2] = size(Xg);
Xs = [reshape(Xg,d1*d2,1),reshape(Yg,d1*d2,1)]';
n = length(Xs);

sigma = 1;
gamma = 3;

kernel = @(x,y)(sigma^2 * exp(-1/2 * gamma *(x - y)'*(x - y)));

R = sqrt(0.5^2 + 0.5^2);
cen = [0.0, 0.0]';
mean = @(x) 0;%1/2/R*((x-cen)'*(x-cen) - R^2);

sigmaNoise = 0.25;
K = [];
for i = 1:m
    for j = 1:m
       K(i,j) = kernel(X(:,i), X(:,j));
       if(i == j)
           K(i,j) = K(i,j) + sigmaNoise*sigmaNoise;
       end
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
       if(i == j)
           Kss(i,j) = Kss(i,j) + sigmaNoise*sigmaNoise;
       end
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
plot(X(1,:), X(2,:), 'r*', 'MarkerSize',6);
surf(Xg,Yg,Fs);
contour(Xg,Yg,Fs,[0 0], 'LineWidth',2,'color', 'r');


figure();
hold on;
plot(X(1,:), X(2,:), 'r*', 'MarkerSize',6);
surf(Xg,Yg,Fs);
contour(Xg,Yg,Fs,[0 0], 'LineWidth',2,'color', 'r');
surf(Xg,Yg,devFsN);
surf(Xg,Yg,devFsP);
axis([-1.5 0 -1.5 1.5])

figure();
hold on;
plot(X(1,:), X(2,:), 'r*');
contour(Xg,Yg,Fs);

% Compute the inside probability.
cdfBell = @(x) 0.5.*(1 + sign(x).*sqrt(1 - exp(-2/pi.*x.*x)));

D = reshape(d, d1,d2);
prob = cdfBell((0-Fs)./D);
figure();
hold on;
surf(Xg,Yg,prob);

