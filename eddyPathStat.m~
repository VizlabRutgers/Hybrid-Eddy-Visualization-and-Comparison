function [outputArg1,outputArg2] = eddyPathStat(eddyHistory,srcData, property,dataFilePath)
% Read data
x_val = double(ncread(srcData, property.x));
y_val = double(ncread(srcData, property.y));
z_val = double(ncread(srcData, property.z));
startLoc = [1,1,1,1];
count = [length(x_val),length(y_val),length(z_val),1];
stride = [1,1,1,1];

startLoc_2D = [1,1,1];
count_2D = [length(x_val),length(y_val),1];
stride_2D = [1,1,1];



% u_val = ncread(srcData, property.u, startLoc, count, stride);
% v_val = ncread(srcData, property.v,startLoc, count, stride);
% w_val = ncread(srcData, property.w,startLoc, count, stride);
% eta_val = ncread(srcData,"ETA",startLoc_2D, count_2D, stride_2D);
% temp_val = ncread(srcData,property.temp,startLoc, count, stride);
% salinity_val = ncread(srcData,"salt",startLoc, count, stride);


for eddyPathIndex = 1:1:size(eddyHistory,1)
    currentPath = eddyHistory{eddyPathIndex};
    
    for eddyFrameIndex = 1:1:size(currentPath,1)
        FrameNum = currentPath(eddyFrameIndex, 15);

        if(currentPath(eddyFrameIndex,14)==1)
            data = load(dataFilePath+"Seperated Structures/clockwise/Frame_"+num2str(FrameNum)+"_eddy_"+num2str(currentPath(eddyFrameIndex,16))+"_statistic.uocd");
        elseif(currentPath(eddyFrameIndex,14)==0)
            data = load(dataFilePath+"Seperated Structures/counterclockwise/Frame_"+num2str(FrameNum)+"_eddy_"+num2str(currentPath(eddyFrameIndex,16))+"_statistic.uocd");
        else
            error('Error: Can not find corresponding eddy data');
        end
        
        currentSurfaceRadius = data(1,13);
        currentVolume = size(data,1);
        currentDepth = currentPath(eddyFrameIndex, 18);
        
        
        
    end  
end
    
    

end

