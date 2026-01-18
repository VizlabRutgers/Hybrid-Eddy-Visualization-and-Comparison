function [eddyData] = readEddyData(dataFilePath, rotationFlag, eddyTimeIndex, objIndex)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    if(rotationFlag==1)
        eddyData = load(dataFilePath+"Seperated Structures/clockwise/Frame_"+num2str(eddyTimeIndex)+"_eddy_"+num2str(objIndex)+"_statistic.uocd");
    elseif(rotationFlag==0)
        eddyData = load(dataFilePath+"Seperated Structures/counterclockwise/Frame_"+num2str(eddyTimeIndex)+"_eddy_"+num2str(objIndex)+"_statistic.uocd");
    else
        error('Error: Can not find corresponding eddy data');
    end
end