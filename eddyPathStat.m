function [eddyPathStat] = eddyPathStat(eddyHistory,srcData, property,dataFilePath)
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



xSize = size(x_val);
ySize = size(y_val);


u_val = ncread(srcData, property.u, startLoc, count, stride);
v_val = ncread(srcData, property.v,startLoc, count, stride);
% w_val = ncread(srcData, property.w,startLoc, count, stride);
eta_val = ncread(srcData,property.eta,startLoc_2D, count_2D, stride_2D);
temp_val = ncread(srcData,property.temp,startLoc, count, stride);
salinity_val = ncread(srcData,property.salt,startLoc, count, stride);

u_val_test = u_val(1:xSize-1,:,:) - u_val(2:xSize,:,:);
u_x = u_val;
u_x(1:xSize-1,:,:) = u_val_test;
u_val_test = u_val(:,1:ySize-1,:) - u_val(:,2:ySize,:);
u_y = u_val;
u_y(:,1:ySize-1,:) = u_val_test;

v_val_test = v_val(1:xSize-1,:,:) - v_val(2:xSize,:,:);
v_x = v_val;
v_x(1:xSize-1,:,:) = v_val_test;
v_val_test = v_val(:,1:ySize-1,:) - v_val(:,2:ySize,:);
v_y = v_val;
v_y(:,1:ySize-1,:) = v_val_test;

s_n = (u_x-v_y).^2;
s_s = (u_y+v_x).^2;
omega = (v_x-u_y).^2;

s_n(s_n==0) = NaN;
s_s(s_n==0) = NaN;
omega(s_n==0) = NaN;
ow_val = s_s+s_n-omega;
temp_val(temp_val==0) = NaN;
salinity_val(salinity_val==0) = NaN;

velocity_mag_val = (sqrt(u_val.^2 +v_val.^2));

eddyPathStat = {};

for eddyPathIndex = 1:1:size(eddyHistory,1)
    currentPath = eddyHistory{eddyPathIndex};
    currentEddyPathStat = cell(size(currentPath,1),16);
    for eddyFrameIndex = 1:1:size(currentPath,1)
        
        FrameNum = currentPath(eddyFrameIndex, 15);

        if(currentPath(eddyFrameIndex,14)==1)
            data = load(dataFilePath+"Seperated Structures/clockwise/Frame_"+num2str(FrameNum)+"_eddy_"+num2str(currentPath(eddyFrameIndex,16))+"_statistic.uocd");
        elseif(currentPath(eddyFrameIndex,14)==0)
            data = load(dataFilePath+"Seperated Structures/counterclockwise/Frame_"+num2str(FrameNum)+"_eddy_"+num2str(currentPath(eddyFrameIndex,16))+"_statistic.uocd");
        else
            error('Error: Can not find corresponding eddy data');
        end
        
        surfaceData = data(data(:,5) == z_val(1),:);
        [~, x_coord_surface] = ismember(round(surfaceData(:,3),4), round(x_val,4));
        [~, y_coord_surface] = ismember(round(surfaceData(:,4),4), round(y_val,4));
        [depthValue, depthValueIndex,~] = unique(data(:,5));
        
        currentSurfaceRadius = data(1,13);
        currentVolume = size(data,1);
        currentDepth = currentPath(eddyFrameIndex, 18);
        currentMinOW = min(data(:,6));
        currentMaxAngular = max(sqrt(data(:,7).^2 + data(:,8).^2)./data(:,13));
        currentMaxSSH = max(abs(eta_val(sub2ind(size(eta_val),x_coord_surface,y_coord_surface))));
        
        currentVorticity = zeros(size(depthValueIndex,1),1);
        currentVelocity = zeros(size(depthValueIndex,1),1);
        currentDepthAverageTempInsideEddy = zeros(size(depthValueIndex,1),1);
        currentDepthBackgroundTemp = zeros(size(depthValueIndex,1),1);
        currentDepthAverageSaltInsideEddy = zeros(size(depthValueIndex,1),1);
        currentDepthBackgroundSalt = zeros(size(depthValueIndex,1),1);
            
        surfaceCenterX = 0;
        surfaceCenterY = 0;
        for depthIndex=1:1:size(depthValueIndex,1)
            centerOnCurrentDepth = data(depthValueIndex(depthIndex),1:2);
            [~, centerX] = ismember(round(centerOnCurrentDepth(1),4), round(x_val,4));
            [~, centerY] = ismember(round(centerOnCurrentDepth(2),4), round(y_val,4));
            [~, centerZ] = ismember(round(data(depthIndex,5),4), round(z_val,4));
            
            if(depthIndex == 1)
                surfaceCenterX = centerOnCurrentDepth(1);
                surfaceCenterY = centerOnCurrentDepth(2);
            end
            currentVorticity(depthIndex) = omega(centerX, centerY, centerZ);
            currentVelocity(depthIndex) = velocity_mag_val(centerX, centerY, centerZ);
            
            currentDepthData = data(data(:,5) == z_val(depthIndex),:);
            currentDepthAverageTempInsideEddy(depthIndex) = mean(currentDepthData(:,11));
            currentDepthBackgroundTemp(depthIndex) = nanmean(nanmean(temp_val(:,:,depthIndex)));
            currentDepthAverageSaltInsideEddy(depthIndex) = mean(currentDepthData(:,10));
            currentDepthBackgroundSalt(depthIndex) = nanmean(nanmean(salinity_val(:,:,depthIndex)));
            
        end
        
        surfaceRadius = unique(surfaceData(:,12));
        currentRadiusVelocity = zeros(length(surfaceRadius)-1,1);
        
        for surfaceRadiusIndex = 2:1:length(surfaceRadius)
            currentRadiusData = surfaceData(surfaceData(:,12)==surfaceRadius(surfaceRadiusIndex),:);
            currentRadiusVelocity(surfaceRadiusIndex-1) = mean(sqrt(currentRadiusData(:,7).^2+currentRadiusData(:,8).^2));
        end
       
        
        currentEddyPathStat(eddyFrameIndex,:) = {currentPath(eddyFrameIndex,15), ...
            surfaceCenterX, surfaceCenterY, ...
            currentDepth, currentSurfaceRadius, currentVolume,  ...
            currentMinOW, currentMaxAngular, currentMaxSSH, ...
            currentVorticity, currentVelocity, ...
            currentDepthAverageTempInsideEddy, currentDepthBackgroundTemp, ...
            currentDepthAverageSaltInsideEddy, currentDepthBackgroundSalt, ...
            currentRadiusVelocity};
        
    end
    
    
    currentEddyPathStat = cell2table(currentEddyPathStat,...
    'VariableNames',{'Frame', 'surfaceCenterX', 'surfaceCenterY',...
    'depth (meter)','surface radius', 'volume',  ...
    'minimum OW', 'maximum angular', 'maximum SSH on surface', ...
    'vorticity along centerline', 'velocity along centerline', ...
    'mean temperature inside eddy on each layer', 'mean bancground temperature on each layer', ...
    'mean salinity inside eddy on each layer', 'mean bancground salinity on each layer', ...
    'mean surface velocity magnitude for all radii(starting by 3)'});

    eddyPathStat{eddyPathIndex}= currentEddyPathStat;
    
end
    
    

end

