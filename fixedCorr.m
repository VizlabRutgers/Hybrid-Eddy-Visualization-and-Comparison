function [corrDistance] = fixedCorr(f1,f2)
%FIXEDCORR 此处显示有关此函数的摘要
%   此处显示详细说明
    if(std(f1) == 0 || std(f2) == 0)
        if(std(f1) == 0 && std(f2) == 0)
            corrDistance = 1;
        else
            corrDistance = 0;
        end
    else
        corrDistance = 1-pdist([f1;f2],'correlation');
    end
end

