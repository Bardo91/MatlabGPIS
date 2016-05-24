% GPIS learning
% 1D example, zero mean

close all; clear all; clc;
% 
% X = [1,0,0;
%     0,1,0;
%     0,0,1;
% %     -0.5774,-0.5774,-0.5774]';
%     -1,-1,-1]';
AppleData;

X = X(:,1:10:end);

X = X(:,X(1,:) < 0);

m = length(X); 

f = zeros(m,1);

[Xg,Yg, Zg] = meshgrid(-0.05:0.01:0.05,-0.05:0.01:0.05,-0.05:0.01:0.05);
[d1,d2] = size(Xg);
Xs = [reshape(Xg,d1*d2,1),reshape(Yg,d1*d2,1),reshape(Zg,d1*d2,1)]';
n = length(Xs);

sigma = 0.75;
gamma = 1;

kernel = @(x,y)(sigma^2 * exp(-1/2 * gamma *(x - y)'*(x - y)));

display('Computing covariance matrix K');
sigmaNoise = 0.01;
K = [];
for i = 1:m
    for j = 1:m
       K(i,j) = kernel(X(:,i), X(:,j));
       if(i == j)
           K(i,j) = K(i,j) + sigmaNoise*sigmaNoise;
       end
    end
end

display('Computing covariance matrix Ks');
Ks = [];
for i = 1:m
    for j = 1:n
       Ks(i,j) = kernel(X(:,i), Xs(:,j)); 
    end
end

display('Computing covariance matrix Kss');
Kss = zeros(n,n);
% 
% if n > 4000
%     delete(gcp)
%     parpool(4)
%     parfor  i = 1:n
%         for j = 1:n
%            Kss(i,j) = kernel(Xs(:,i), Xs(:,j));
%            if(i == j)
%                Kss(i,j) = Kss(i,j) + sigmaNoise*sigmaNoise;
%            end
%         end
%     end
% else
   for  i = 1:n
         for j = 1:n
           Kss(i,j) = kernel(Xs(:,i), Xs(:,j));
           if(i == j)
               Kss(i,j) = Kss(i,j) + sigmaNoise*sigmaNoise;
           end
        end
    end 

% end

figure();
imagesc(K);
colorbar;

figure();
imagesc(Ks);
colorbar;

figure();
imagesc(Kss);
colorbar;

R = 0.3;
cen = [0.0, 0.0, 0.0]';
mean = @(x) 1/2/R*((x-cen)'*(x-cen) - R^2);

display('Computing mean vectors');
mu = zeros(m,1);
for i = 1:m
    mu(i) = mean(X(:,i));
end

mus = zeros(n,1);
for i = 1:n
    mus(i) = mean(Xs(:,i));
end

display('Computing regression');
kinv = inv(K);
fs = mus + Ks'*kinv*(f - mu);
sig = Kss' - Ks'*kinv*Ks;

figure();
imagesc(sig);
colorbar;

Fs = reshape(fs,d1,d1,d1);
figure();
hold on;
plot3(X(1,:), X(2,:), X(3,:), 'r.', 'MarkerSize',40);
p = patch(isosurface(Xg, Yg, Zg, Fs, 0));
p.FaceColor = 'green';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3)
camlight; lighting phong;


d = diag(sig);
devFsP = fs - d;
devFsP = reshape(devFsP,d1,d2);

devFsN = fs + d;
devFsN = reshape(devFsN,d1,d2);

devFsN = reshape(devFsN,d1,d1,d1);
devFsP = reshape(devFsP,d1,d1,d1);

figure();
hold on;
pMean = patch(  isosurface(Xg, Yg, Zg, Fs, 0), ...
                'FaceColor','green',...
                'FaceAlpha',0.5,...
                'EdgeColor', 'none');

pDevP = patch(  isosurface(Xg, Yg, Zg, devFsP, 0), ...
                'FaceColor','red',...
                'FaceAlpha',0.5,...
                'EdgeColor', 'none');

pDevN = patch(  isosurface(Xg, Yg, Zg, devFsN, 0), ...
                'FaceColor','blue',...
                'FaceAlpha',0.5,...
                'EdgeColor', 'none');
            
daspect([1 1 1])
view(3)
plot3(X(1,:), X(2,:), X(3,:), 'r.', 'MarkerSize',40);
camlight; lighting phong;


display('Computing the inside probability');
% Compute the inside probability.
cdfBell = @(x) 0.5.*(1 + sign(x).*sqrt(1 - exp(-2/pi.*x.*x)));

D = reshape(d, d1,d1,d1);
prob = cdfBell((0-Fs)./D);
figure();
hold on;
p0 = patch(  isosurface(Xg, Yg, Zg, prob, 0), ...
                'FaceColor','cyan',...
                'FaceAlpha',0.25,...
                'EdgeColor', 'none');
p25 = patch(  isosurface(Xg, Yg, Zg, prob, 0.25), ...
                'FaceColor','blue',...
                'FaceAlpha',0.25,...
                'EdgeColor', 'none');
p50 = patch(  isosurface(Xg, Yg, Zg, prob, 0.5), ...
                'FaceColor','green',...
                'FaceAlpha',0.25,...
                'EdgeColor', 'none');
p75 = patch(  isosurface(Xg, Yg, Zg, prob, 0.75), ...
                'FaceColor','red',...
                'FaceAlpha',0.25,...
                'EdgeColor', 'none');
p1 = patch(  isosurface(Xg, Yg, Zg, prob, 1), ...
                'FaceColor','yellow',...
                'FaceAlpha',0.5,...
                'EdgeColor', 'none');
daspect([1 1 1])
plot3(X(1,:), X(2,:), X(3,:), 'r.', 'MarkerSize',40);
view(3)
camlight; lighting phong;

figure()
hold on;
plot3(X(1,:), X(2,:), X(3,:), 'r.', 'MarkerSize',40);
p = patch(isosurface(Xg, Yg, Zg, Fs, 0),'FaceAlpha',0.7);
p.FaceColor = 'green';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3)
camlight; lighting phong;
colormap('default');   % set colormap
contourf(Xg(:,:,1),Yg(:,:,1),prob(:,:,floor(d1/2)));        % draw image and scale colormap to values range
colorbar;