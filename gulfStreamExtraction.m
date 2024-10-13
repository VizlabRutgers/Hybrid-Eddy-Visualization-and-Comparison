function [streamRidge_result] = gulfStreamExtraction(srcData,property,spaceLimit)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
x_val = double(ncread(srcData, property.x));
y_val = double(ncread(srcData, property.y));
z_val = double(ncread(srcData, property.z));
time_val = double(ncread(srcData, property.time));

[~, lower_x] = min(abs(x_val - spaceLimit.x0));
[~, upper_x] = min(abs(x_val - spaceLimit.x1));
[~, lower_y] = min(abs(y_val - spaceLimit.y0));
[~, upper_y] = min(abs(y_val - spaceLimit.y1));

streamRidge_result = zeros(upper_x-lower_x+1,length(time_val));

for FrameNum = 1:1:length(time_val)


    
    startLoc_2D = [lower_x,lower_y,1,FrameNum];
    count_2D = [upper_x-lower_x+1,upper_y-lower_y+1,1,1];
    stride_2D = [1,1,1,1];
    
    
    
    u_val = ncread(srcData, property.u, startLoc_2D, count_2D, stride_2D);
    v_val = ncread(srcData, property.v,startLoc_2D, count_2D, stride_2D);
    
    vel_mag = sqrt(u_val.^2+v_val.^2);
    vel_mag(isnan(vel_mag)) = 0;
    
    
    % streamRidge = im2uint8(watershed(vel_mag,8));
    
    

    
    
    vel_mag_cropped = vel_mag;
    
    
    x_val_cropped = x_val(lower_x:upper_x);
    y_val_cropped = y_val(lower_y:upper_y);
    
    vel_mag_hmin = imhmin(vel_mag_cropped,0.3);
    
    streamRidge = watershed(vel_mag_hmin,8);
    
    streamRidge = ~(imbinarize(streamRidge,0));
    streamRidge(vel_mag_cropped==0) = false;



    [ridgeRow, ridgeCol] = find(streamRidge==1);

    ridge_XY = zeros(length(x_val_cropped),1);

    for xIndex = 1:1:length(ridgeRow)
        ridgeX = ridgeRow(ridgeRow(xIndex));
        ridgeY = ridgeCol(xIndex);

        if(ridge_XY(ridgeRow(xIndex))==0)
            ridge_XY(ridgeRow(xIndex)) = ridgeY;
        else
            if(ridge_XY(ridgeRow(xIndex))>ridgeY)
                ridge_XY(ridgeRow(xIndex)) = ridgeY;
            end
        end
    end

    streamRidge_result(:,FrameNum) = ridge_XY;
    
    figure,
    uimagesc(gca,x_val_cropped, y_val_cropped, vel_mag_cropped');
    hold on;
    set(gca, 'YDir', 'normal');
    % daspect([1 1 1])
    

    plot(x_val_cropped(ridgeRow), y_val_cropped(ridgeCol), 'r.', 'MarkerSize',10);

end
end