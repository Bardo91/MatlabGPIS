close all
clear variables
clc;

% load the part location and surface normal data (thsi would be the output from an oversegmentation) 
% load  '/Users/Yannick/Google Drive/PDM/OptHyperparam/simulated_data/Simple side high resolution/bmw_11'
% load '/Users/Yannick/Google Drive/PDM/OptHyperparam/Ellipse_plot/car_simul_cleandata'
load '/Users/Yannick/Google Drive/PDM/OptHyperparam/bmw_total'

% initialise the hyperparameters (always needed for iterative method)
noiseVals0 = .001;
noiseGrad0 = .037;
sigma0 = 0.16;%0.2844; %0.1
L0 = .2823; 
gamma0 = 2.24;%0.6525;%1/L0^2;
cx0 =    -0.1075;%-0.12;%-6;
cy0 =    0.1002;%0.227;%0;
cz0 =    -1.7149;%-1.34;%-3.1;
a0 = 1.0721;%0.98;%0.97;
b0 = 2.2142;%2.95;%2.65;
c0 = 1.4802;%1.21;%1.29;

% Optimization variable 1x3 vector. 
% - Optim_var(1) if optim sigma and gamma.
% - Optim_var(2) if optim cx cy cz.
% - Optim_var(3) if optim a,b,c.
Optim_var = [1 1 1];



[sigma,gamma,noiseVals,noiseGrad,cx,cy,cz,a,b,c] = ...
    FindHyperparamRed(PartMeans,SurfNormals,sigma0,gamma0,noiseVals0,noiseGrad0,cx0,cy0,cz0,a0,b0,c0,Optim_var)


%% Value of optimization full complete car

% sigma0 = 0.16;%0.2844; %0.1
% L0 = .2823; 
% gamma0 = 2.24;%0.6525;%1/L0^2;
% cx0 =    -0.1075;%-0.12;%-6;
% cy0 =    0.1002;%0.227;%0;
% cz0 =    -1.7149;%-1.34;%-3.1;
% a0 = 1.0721;%0.98;%0.97;
% b0 = 2.1909;%2.95;%2.65;
% c0 = 1.4733;%1.21;%1.29;



%% These values are for the bmw_11/25/...etc
% cx0 =    -6.2643;%-6;
% cy0 =    0.1404;%0;
% cz0 =    -3.3735;%-3.1;
% a0 = 1.38;
% b0 = 2.44;
% c0 = 1.28;