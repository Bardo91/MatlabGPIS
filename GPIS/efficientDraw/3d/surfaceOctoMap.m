function [ surfacePoints ] = surfaceOctoMap( prior, dataPoints, preExpandingIterations, iterations, debugMode)
%surfaceOctoMap Summary of this function goes here
%   Detailed explanation goes here

    X = dataPoints(1:3, :);
    m = length(X); 
    
    f = zeros(1,length(dataPoints));
    f = [f ; dataPoints(4:6,:)]; 
    f = reshape(f, [1,4*length(X)])';

    sigma = prior.sigma;
    gamma = prior.gamma;
    noiseVal = prior.noiseVal;
    noiseGrad = prior.noiseGrad;
    mean = prior.mean;
    
    if(debugMode)
        display('Computing covariances')
    end
    K = ComputeFullKder(sigma, gamma, X, noiseVal, noiseGrad);
    
    
    if(debugMode)
        display('Computing means')
    end
    for i = 1:m
        mu(((i-1)*4+1):((i-1)*4 +4)) = mean(X(:,i));
    end

    %% Efficiend draw
    xLimits = [min(X(1,:)), max(X(1,:))];
    xLimits(1) = xLimits(1)*1.2;
    xLimits(2) = xLimits(2) + abs(xLimits(2))*1.2;
    yLimits = [min(X(2,:)), max(X(2,:))];
    yLimits(1) = yLimits(1)*1.2;
    yLimits(2) = yLimits(2) + abs(yLimits(2))*1.2;
    zLimits = [min(X(3,:)), max(X(3,:))];
    zLimits(1) = zLimits(1)*1.2;
    zLimits(2) = zLimits(2) + abs(zLimits(2))*1.2;
    centroid = [sum(xLimits)/2, sum(yLimits)/2,sum(zLimits)/2];

    % 2 valid, 1 new, 0 invalid.
    root = {2,  {}, xLimits, yLimits, centroid, -1, zLimits};
    Qmat = K\(f - mu');

    evalFun = @(x) mean(x) + ComputeKderX1X2(sigma, gamma, x, X)*Qmat;

    
    if(debugMode)
        display('Computing map')
    end
    
    
    if(debugMode)
        display('Preexpanding cells')
    end
    
    for i =1:preExpandingIterations
        if(debugMode)
            display(['Preexpanding iteration ', num2str(i)])
        end
        root = expandCell(root,evalFun, 2);
    end
    
    for iter=1:iterations
        if(debugMode)
            display(['Computing iteration ', num2str(iter)])
        end
        % Expand tree
        root = expandCell(root, evalFun, 1); 
        % Validate branches
        root = validatePoints(root, root, evalFun);
    end
    root = validatePoints(root, root, evalFun);
    surfacePoints  = [];
    cols = [];
    [surfacePoints, cols] = getPointsTree(surfacePoints, cols, root, true);
end

