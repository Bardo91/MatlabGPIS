function [ probs ] = ProbabilityMap( prior, dataPoints, dataNormals, evalPoints )
    %PROBABILITYMAP Summary of this function goes here
    %   Detailed explanation goes here

    nData = length(dataPoints);
    nEval = length(evalPoints);

    dims = length(dataPoints(:,1));

    %% Compute means
    mu = zeros(nData*(dims+1));
    for i = 1:nData
        mu((i-1)*4 +1) = mean(dataPoints(:,i));
        mu((i-1)*4 +2) = meandx(dataPoints(:,i));
        mu((i-1)*4 +3) = meandy(dataPoints(:,i));
        mu((i-1)*4 +4) = meandz(dataPoints(:,i));
    end

    mus = zeros(nEval*(dims+1));
    for i = 1:nEval
        is = (i-1)*(dims+1);
        ie = (i-1)*(dims+1) + (dims+1);
        mus(is:ie) = prior.mean(evalPoints(:,i));
    end


    %% Compute Kernels
    K = zeros(nData*(dims+1), nData*(dims+1));
    for i=1:nData
       for j=1:nData
            K(i,j) = prior.kernel(dataPoints(:,i), dataPoints(:,j));
       end
    end

    Ks = zeros(nData*(dims+1), nEval*(dims+1));
    for i=1:nData
       for j=1:nEval
            K(i,j) = prior.kernel(dataPoints(:,i), evalPoints(:,j));
       end
    end

    Kss = zeros(nEval*(dims+1))
    
    %% Regression
    kinv = inv(K);
    fs = mus + Ks*kinv*(f - mu);  
    sig = Kss - diag(Ks*kinv*Ks');

    %% Probabilities.
    cdfBell = @(x) 0.5.*(1 + sign(x).*sqrt(1 - exp(-2/pi.*x.*x)));

    probs = cdfBell((0-fs)./d);

end

