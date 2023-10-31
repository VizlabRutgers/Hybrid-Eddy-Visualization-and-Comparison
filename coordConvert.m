function [coordOutput] = coordConvert(coordInput, data_val, mode, resolution)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if(mode == 1)
    coordOutput = data_val(coordInput +1);
elseif(mode == 2)
    [~, coordOutput] = ismember(round(coordInput,4),round(data_val,4));
elseif(mode == 3)
    coordOutput = double((coordInput-data_val(1))./resolution+1);
end