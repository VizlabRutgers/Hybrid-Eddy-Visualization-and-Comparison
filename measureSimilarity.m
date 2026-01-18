function [euclideanDistance,cosineDistance,corrDistance,ssimDistance] = measureSimilarity(f1, f2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

length_1 = length(f1);
length_2 = length(f2);

minLen = min(length_1, length_2);

if(length(f1) == length(f2))
    euclideanDistance = pdist([f1;f2]);
    cosineDistance = 1-pdist([f1;f2],'cosine');
    if(minLen>3)
        corrDistance = movcorr(f1,f2,3);  
    else
        corrDistance = fixedCorr(f1,f2);
    end
    ssimDistance = ssim(f1',f2');
elseif(length(f1) > length(f2))
    f2_fillZero = [f2, zeros(1,length(f1) - length(f2))];
    scaleRatio = length(f1)/length(f2);
    f1_scaled_interp = interp1((0:1:length(f1)-1)/scaleRatio, f1/scaleRatio, 0:1:(length(f2)-1));

    euclideanDistance = pdist([f1;f2_fillZero]);
    cosineDistance = 1-pdist([f1_scaled_interp;f2],'cosine');

    if(minLen>3)
        corrDistance = movcorr(f1,f2_fillZero,3);  
    else
        corrDistance = fixedCorr(f1,f2_fillZero);
    end

    ssimDistance = ssim(f1,f2_fillZero);
else
    f1_fillZero = [f1, zeros(1,length(f2) - length(f1))];
    scaleRatio = length(f2)/length(f1);
    f2_scaled_interp = interp1((0:1:length(f2)-1)/scaleRatio, f2/scaleRatio, 0:1:(length(f1)-1));

    euclideanDistance = pdist([f1_fillZero;f2]);
    cosineDistance = 1-pdist([f1;f2_scaled_interp],'cosine');

    if(minLen>3)
        corrDistance = movcorr(f1_fillZero,f2,3);  
    else
        corrDistance = fixedCorr(f1_fillZero,f2);
    end

    ssimDistance = ssim(f1_fillZero,f2);

end
end