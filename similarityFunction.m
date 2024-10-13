function [euclideanDistance,cosineDistance,corrDistance,ssimDistance] = similarityFunction(f1,f2)
%SIMILARITYFUNCTION Summary of this function goes here
%   Detailed explanation goes here
euclideanDistance = pdist([f1;f2]);
cosineDistance = 1-pdist([f1;f2],'cosine');
corrDistance = 1-pdist([f1;f2],'correlation');
ssimDistance = ssim(f1,f2);
end