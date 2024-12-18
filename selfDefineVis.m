function [outputArg1,outputArg2] = selfDefineVis(allEddy, frameTotal, eddyPathIndex, srcData, property,dataFilePath, eddyPathHistory)
%SELFDEFINEVIS Summary of this function goes here
%   Detailed explanation goes here
    close all;
    x_val = double(ncread(srcData, property.x));
    y_val = double(ncread(srcData, property.y));
    z_val = double(ncread(srcData, property.z));
    individualEddyHistory = eddyPathHistory{eddyPathIndex};
    
    for frameIndex = frameTotal
    
    startLoc = [1,1,1,frameIndex];
    count = [length(x_val),length(y_val),length(z_val),1];
    stride = [1,1,1,1];
%     startLoc_2D = [1,1,FrameNum];
%     count_2D = [length(x_val),length(y_val),1];
%     stride_2D = [1,1,1];

    u_val = ncread(srcData, property.u, startLoc, count, stride);
    v_val = ncread(srcData, property.v,startLoc, count, stride);
%     w_val = ncread(srcData, "W",startLoc, count, stride);
%     eta_val = ncread(srcData,"ETA",startLoc_2D, count_2D, stride_2D);
    temp_val = ncread(srcData,property.temp,startLoc, count, stride);
    
    

    
        individualEddy = individualEddyHistory(frameIndex,:);

        if(individualEddyHistory(frameIndex,14)==1)
            data = load(dataFilePath+"Seperated Structures/clockwise/Frame_"+num2str(individualEddyHistory(frameIndex,15))+"_eddy_"+num2str(individualEddyHistory(frameIndex,16))+"_statistic.uocd");
        elseif(individualEddyHistory(frameIndex,14)==0)
            data = load(dataFilePath+"Seperated Structures/counterclockwise/Frame_"+num2str(individualEddyHistory(frameIndex,15))+"_eddy_"+num2str(individualEddyHistory(frameIndex,16))+"_statistic.uocd");
        else
            error('Error: Can not find corresponding eddy data');
        end

        x = data(:,3);
        y = data(:,4);
        z = data(:,5);
        ow = data(:,6);
        u_obj = data(:,7);
        v_obj = data(:,8);
        temp = data(:,10);
        salt = data(:,11);
        
        maxZ = max(z(:));
        
        [~,maxZ_coord] = ismember(round(maxZ,2), round(z_val,2));

        bound = boundary(x,y,z,0.9);
        figure,
        localOb = trimesh(bound,x,y,z,sqrt(u_obj.^2+v_obj.^2), 'FaceAlpha', '0.5', 'EdgeAlpha', 0.8);
        set(gca, "ZDir", 'reverse');
        
        title("eddy in frame " + num2str(frameIndex));
        
        figure,
        quiver(x_val, y_val, u_val(:,:,1)', v_val(:,:,1)');
        daspect([1,1,1])
        
        title("surface velocity field in frame " + num2str(frameIndex));
        
        
        figure,
        quiver(x_val, y_val, u_val(:,:,maxZ_coord+1)', v_val(:,:,maxZ_coord+1)', 2);
        daspect([1,1,1])
        
        title("bottom velocity field in frame " + num2str(frameIndex));
    end

end

