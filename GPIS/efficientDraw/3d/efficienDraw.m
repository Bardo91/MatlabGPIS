% GPIS learning
% 3D example, Sphere mean and normals

close all; clear all; clc;

load('bmw_11.mat');

data = [PartMeans;SurfNormals];

prior.sigma = 0.12;
prior.gamma = 2.6;
prior.noiseVal = 0.001;
prior.noiseGrad = 0.001;

R = 1;
cen = (sum(PartMeans')/length(PartMeans))';
prior.mean = @(x) [   1/2/R*((x-cen)'*(x-cen) - R^2);...
                1/R*((x(1)-cen(1)));...
                1/R*((x(2)-cen(2)));...
                1/R*((x(3)-cen(3)))];
            
surface = surfaceOctoMap(prior, data, 5, true);

figure();
hold on;
points = [];
cols = [];
plot3(surface(:,1), surface(:,2),surface(:,3), 'go');
plot3(data(1,:), data(2,:), data(3,:), 'r.', 'MarkerSize',40);
quiver3(data(1,:), data(2,:), data(3,:), data(4,:), data(5,:),data(6,:));
grid;
axis equal


% Ground Truth
% [Xg,Yg, Zg] = meshgrid(-2:0.5:2,-1.4:0.5:2,-2:0.5:2);
% [d1,d2] = size(Xg);
% Xs = [reshape(Xg,d1*d2,1),reshape(Yg,d1*d2,1),reshape(Zg,d1*d2,1)]';
% n = length(Xs);
% 
% for i = 1:n
%     mus((i-1)*4 +1) = mean(Xs(:,i));
%     mus((i-1)*4 +2) = meandx(Xs(:,i));
%     mus((i-1)*4 +3) = meandy(Xs(:,i));
%     mus((i-1)*4 +4) = meandz(Xs(:,i));
% end
% mus = mus';
% Ks = ComputeKderX1X2(sigma, gamma, Xs, X)';
% Kss = ComputeFullKder(sigma, gamma, Xs, noiseVal, noiseGrad);
% 
% fs = mus + Ks'*inv(K)*(f - mu');
% sig = Kss' - Ks'*inv(K)*Ks;
% 
% Fs = reshape(fs(1:4:end),7,9,9);
% 
% figure();
% hold on;
% plot(X(1,:), X(2,:), 'r.', 'MarkerSize',40);
% figure();
% hold on;
% pMean = patch(  isosurface(Xg, Yg, Zg, Fs, 0), ...
%                 'FaceColor','green',...
%                 'FaceAlpha',0.5,...
%                 'EdgeColor', 'none');
% plot3(X(1,:), X(2,:), X(3,:), 'r.', 'MarkerSize',40);
% quiver3(X(1,:), X(2,:), X(3,:), data(2,:), data(3,:),data(4,:));