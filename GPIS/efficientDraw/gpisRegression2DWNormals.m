% GPIS learning
% 2D example, Sphere mean and normals

close all; clear all; clc;

X = [0,-0.5;
    -0.3,-0.1;
    0.5,-1.0;
     -0.5,0.5;
     0.5,0.5]';

f = [   0,-cos(20/180*pi),-sin(20/180*pi),...
        0,-cos(45/180*pi),-sin(45/180*pi),...
        0,-cos(120/180*pi),-sin(120/180*pi),...
        0, -cos(45/180*pi), sin(45/180*pi),...
        0, cos(45/180*pi), sin(45/180*pi)]';

data = reshape(f, [3,length(X)])';
m = length(X); 

sigma = 1;
gamma = 1;

display('Computing means');
R = sqrt(0.5^2 + 0.5^2);
cen = [0.0, 0.0]';
mean = @(x) 1/2/R*((x-cen)'*(x-cen) - R^2);
meandx = @(x) 1/2/R*(2*(x(1)-cen(1)));
meandy = @(x) 1/2/R*(2*(x(2)-cen(2)));

display('Computing mean and covariance of data');
K = ComputeFullKder(sigma, gamma, X, 0.2, 0);
for i = 1:m
    mu((i-1)*3 +1) = mean(X(:,i));
    mu((i-1)*3 +2) = meandx(X(:,i));
    mu((i-1)*3 +3) = meandy(X(:,i));
end
mu = mu';

iterations = 6;

xLimits = [-2.0, 2.0];
yLimits = [-2.0, 2.0];
centroid = [sum(xLimits)/2, sum(yLimits)/2];

% 2 valid, 1 new, 0 invalid.
root = {2,  {}, xLimits, yLimits, centroid, -1};
Qmat = inv(K)*(f - mu);

evalFun = @(x) [mean(x), meandx(x), meandy(x)]' + ComputeKderX1X2(sigma, gamma, x, X)*Qmat;

for iter=1:iterations
    % Expand tree
    root = expandCell(root,evalFun);
    
    % Validate branches
    
end

% [Xg,Yg] = meshgrid(-1.4:0.2:1.4,-1.4:0.2:1.4);
% [d1,d2] = size(Xg);
% Xs = [reshape(Xg,d1*d2,1),reshape(Yg,d1*d2,1)]';
% n = length(Xs);

% display('Computing covariance matrix Ks');
% Ks = ComputeKderX1X2(sigma, gamma, Xs, X);
% Ks = Ks';
% 
% for i = 1:n
%     mus((i-1)*3 +1) = mean(Xs(:,i));
%     mus((i-1)*3 +2) = meandx(Xs(:,i));
%     mus((i-1)*3 +3) = meandy(Xs(:,i));
% end
% 
% mus = mus';
% 
% display('Computing Regression');
% fs = mus + Ks'*inv(K)*(f - mu);
% sig = Kss' - Ks'*inv(K)*Ks;

display('displaying')
figure();
hold on;
plot(X(1,:), X(2,:), 'r.', 'MarkerSize',40);
points = [];
cols = [];
[points, cols] = getPointsTree(points, cols, root);
plot3(points(:,1), points(:,2), cols, 'o')
quiver(X(1,:)', X(2,:)', data(:,2), data(:,3));
grid;
axis([-2 2 -2 2]);
