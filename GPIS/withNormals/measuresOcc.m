function [diffSqrt, corrOcc] = measuresOcc(groundTruth, map) 
    [n1, m1] = size(groundTruth);
    [n2, m2] = size(map);

    A = reshape(groundTruth,[1,n1*m1]);
    B = reshape(map,[1,n2*m2]);
    
    diffSqrt = sum(sum((A-B).^2));
    corrOcc = (mean(A.*B) - mean(A)*mean(B))/(sqrt(var(A))*sqrt(var(B)));
end