function [propertyEvaluationResult] = eddyCoordEvaluationFunc(interpProperty,groundtruthProperty)
%ELLIPSEEVALUATIONFUNC Summary of this function goes here
%   Detailed explanation goes here
dtwDistance = dtw(groundtruthProperty,interpProperty);
    
[euclideanDistance,cosineDistance,corrDistance,ssimDistance] = measureSimilarity(interpProperty', groundtruthProperty');
propertyEvaluationResult = [euclideanDistance,cosineDistance,corrDistance,ssimDistance,dtwDistance];

end

