X = PartMeans;
norms = SurfNormals;
% norms = norms(:,X(1,:) < 0);
% X = X(:,X(1,:) < 0);
m = length(X);

f = zeros(m,1);
f = [f,norms'];
[a b] = size(f);
f = reshape(f', a*b, 1);

step = 0.1;
limx=3.5;
limy = 3.2;
%  ellipsoid mean
sigma = Prior.Sigma; 
gamma = Prior.Gamma; 
noiseVals = Prior.noiseVals; 
noiseGrad = Prior.noiseGrad; 
prior = Prior;       

[meanValue, meanGrad] = computePriorFunctions(Prior);

% %  spherical mean
% sigma = 0.4; 
% gamma = 10; 
% noiseVals = 0.06; 
% noiseGrad = 0.04; 
% R = 1;
% cen = [0.0, 0.0, 0.0]';
% mean = @(x)1/2/R*((x-cen)'*(x-cen) - R^2);
% meandx = @(x) 1/R*((x(1)-cen(1)));
% meandy = @(x) 1/R*((x(2)-cen(2)));
% meandz = @(x) 1/R*((x(3)-cen(3)));
% prior = struct( 'pos',[0 0 0],'type','S',...
%                 'param', [R R R], 'rot', [2*pi 0 0],...
%                 'Sigma', sigma, 'Gamma', gamma,...
%                 'noiseVals', noiseVals, 'noiseGrad', noiseGrad);
%          
% % % No prior
% sigma = 0.11;%0.5395;
% gamma = 11.2;%10.43;    
% noiseVals = 1.17e-5;
% noiseGrad = 5.87e-4;
% meanLevel = 0.17;
% mean = @(x)meanLevel;
% meandx = @(x)0;
% meandy = @(x)0;
% meandz = @(x)0;
% prior = struct( 'pos',[0 0 0],'type','N',...
%                 'param', [meanLevel 1 1], 'rot', [2*pi 0 0],...
%                 'Sigma', sigma, 'Gamma', gamma,...
%                 'noiseVals', noiseVals, 'noiseGrad', noiseGrad);
% 
% 

[Xg,Yg,Zg] = meshgrid(-limx:step:limx,-limy:step:limy,[0]);
[d1,d2,d3] = size(Xg);
Xs = [reshape(Xg,d1*d2*d3,1),reshape(Yg,d1*d2*d3,1),reshape(Zg,d1*d2*d3,1)]';
n = length(Xs);

display('Computing covariance matrix K');
K = ComputeFullKder(sigma, gamma, X, noiseVals, noiseGrad);

display('Computing covariance matrix Ks');
Ks = ComputeKderX1X2(sigma, gamma, X, Xs);

display('Computing covariance matrix Kss');
[D,N] = size(Xs);
Kss = DiagComputeKderX1X2(sigma,gamma,Xs,Xs);
NoiseDiag =  [noiseVals * ones(1,N);
              noiseGrad * ones(D,N)];
Kss = Kss + diag(NoiseDiag(:));
% Kss = ComputeFullKder(sigma, gamma, Xs, noiseVals, noiseGrad);

display('Computing means');

mu = zeros(m*4,1);
for i = 1:m
    mu((i-1)*4 +1) = meanValue(X(:,i));
    mu((i-1)*4 +2:(i-1)*4 +4) = meanGrad(X(:,i));
end
mus = zeros(n*4,1);
for i = 1:n
    mus((i-1)*4 +1) = meanValue(Xs(:,i));
    mus((i-1)*4 +2:(i-1)*4 +4) = meanGrad(Xs(:,i));
end

display('Computing regression');
R = K\Ks;
fs = mus + R'*(f - mu);
var = diag(Kss) - diag(R'*Ks);

sig = sqrt(var);
SIG = reshape(sig(1:4:end), d1,d2,d3);
Fs = reshape(fs(1:4:end),d1,d2,d3);

display('Computing the inside probability');
% Compute the inside probability.
cdfBell = @(x) 0.5*(1 + sign(x).*sqrt(1 - exp(-2/pi*x.*x)));
% prob = cdfBell((0-Fs)./SIG);
prob = normcdf(0, Fs, SIG);
%% Accurate plot
[faces, vertices] = computeSurface(X, norms, prior, meanValue, meanGrad, X(:,1), 0.2, false);

figure
hold on
axis equal
vertices(:,3) = - vertices(:,3);
quiver3(X(1,:),X(2,:),X(3,:), norms(1,:),norms(2,:),norms(3,:),'linewidth',2,'color','r');

patch('faces',faces,'vertices',vertices,...
    'facecolor',[0.5 0.5 0.5], ...
    'edgecolor', 'none', ...
    'facelighting','phong',...
    'FaceAlpha', 1);
camlight
set(gca,'view',[46.8000   18.8000]);
light('Position',[-1 -1 0])
view([25 30])
V = [   -limx,-limy,0;
        -limx,limy,0;
        limx,limy,0;
        limx,-limy,0];
F = [1,2,3;
    1,3,4]
patch('Vertices',V,'Faces',F,...
    'facecolor',[0.5 0.3 0.3], ...
    'edgecolor', 'none', ...
    'facelighting','phong',...
    'FaceAlpha', 0.3)

figure
hold on
axis equal
quiver3(X(1,:),X(2,:),X(3,:), norms(1,:),norms(2,:),norms(3,:),'linewidth',2,'color','r');
patch('faces',faces,'vertices',vertices,...
    'facecolor',[0.5 0.5 0.5], ...
    'edgecolor', 'none', ...
    'facelighting','phong',...
    'FaceAlpha', 0.5);
camlight
set(gca,'view',[46.8000   18.8000]);
light('Position',[-1 -1 0])
contourf(Xg(:,:,1),Yg(:,:,1),prob(:,:,1));        % draw image and scale colormap to values range
colorbar;
caxis([0 1])
view([25 30])

figure
% subplot(1,4,2)
imagesc(prob(:,:,1));        % draw image and scale colormap to values range
colorbar;
caxis([0 1])
axis equal
% 
% figure
% % subplot(1,4,3)
% contourf(Xg(:,:,1),Yg(:,:,1),prob(:,:,2));        % draw image and scale colormap to values range
% colorbar;
% axis equal
% 
% figure
% % subplot(1,4,4)
% contourf(Xg(:,:,1),Yg(:,:,1),prob(:,:,3));        % draw image and scale colormap to values range
% colorbar;
% axis equal
