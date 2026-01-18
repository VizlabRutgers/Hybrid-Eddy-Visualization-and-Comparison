function [outputArg1,outputArg2] = selfDefineVis(allEddy, frameTotal, eddyPathIndex, srcData, property,dataFilePath, eddyPathHistory)
%SELFDEFINEVIS Summary of this function goes here
%   Detailed explanation goes here
    close all;
    x_val = double(ncread(srcData, property.x));
    y_val = double(ncread(srcData, property.y));
    z_val = double(ncread(srcData, property.z));
    individualEddyHistory = eddyPathHistory{eddyPathIndex};
    
    minDepth = inf;
    minDepthTemp = 0;
    
    startFrame = individualEddyHistory(1,15);
    for frameIndex = frameTotal
        

        
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
        
        minX = min(x);
        maxX = max(x);
        minY = min(y);
        maxY = max(y);
        
        maxZ = max(z(:));
        
        centerX = data(1,1);
        centerY = data(1,2);
        
        [~,centerX_coord] = ismember(round(centerX,2), round(x_val,2));
        [~,centerY_coord] = ismember(round(centerY,2), round(y_val,2));
        
        startLoc = [1,1,1,frameIndex+startFrame-1];
        count = [length(x_val),length(y_val),length(z_val),1];
        stride = [1,1,1,1];
    %     startLoc_2D = [1,1,FrameNum];
    %     count_2D = [length(x_val),length(y_val),1];
    %     stride_2D = [1,1,1];

        u_val = ncread(srcData, property.u, startLoc, count, stride);
        v_val = ncread(srcData, property.v,startLoc, count, stride);
    %     w_val = ncread(srcData, "W",startLoc, count, stride);
    %     eta_val = ncread(srcData,"ETA",startLoc_2D, count_2D, stride_2D);
%         temp_val = ncread(srcData,property.temp,startLoc, count, stride);


        velocity_mag_val = sqrt(u_val.^2 + v_val.^2);

        xSize = size(x_val);
        ySize = size(y_val);
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
%         ow_val = ow_val./(velocity_mag_val.^2);
%         colorRange = [-0.05,0];
        colorRange = [-0.0005,0];



        [~,maxZ_coord] = ismember(round(maxZ,2), round(z_val,2));
        minDepthTemp = maxZ_coord;
        if(minDepth>minDepthTemp)
            minDepth = minDepthTemp;       
        end
        

        maxZ_coord = minDepth+1;
        
        bound = boundary(x,y,z,0.9);
        figure,
        localOb = trimesh(bound,x,y,z,sqrt(u_obj.^2+v_obj.^2), 'FaceAlpha', '0.5', 'EdgeAlpha', 0.8);
        set(gca, "ZDir", 'reverse');

        title("eddy in frame " + num2str(frameIndex));

%         figure,
%         quiver(x_val, y_val, u_val(:,:,1)', v_val(:,:,1)');
%         daspect([1,1,1])
% 
%         title("surface velocity field in frame " + num2str(frameIndex));


        figure, imagesc(y_val(1:70), z_val, squeeze(ow_val(centerX_coord,1:70,:))');
        colorbar();
        caxis([-0.003,0])
        title("vertical OW clip in frame " + num2str(frameIndex));
        xlabel("latitude");
        ylabel("depth");
        
        figure, imagesc(y_val(1:70), z_val, squeeze(velocity_mag_val(centerX_coord,1:70,:))');
        colorbar();
        caxis([0,0.4])
        title("vertical Velocity Magnitude clip in frame " + num2str(frameIndex));
        xlabel("latitude");
        ylabel("depth");
        
        
        figure,
        uimagesc(gca,x_val, y_val, ow_val(:,:,maxZ_coord-1)');
        cb = colorbar();
        caxis(colorRange);
        hold on

        quiver(x_val, y_val, u_val(:,:,maxZ_coord-1)', v_val(:,:,maxZ_coord-1)', 5,'r');
        daspect([1,1,1])
        hold on
        title("bottom - 1 layer velocity field in frame " + num2str(frameIndex));

        

        dataInThisDepth = data(round(data(:,5),1) == round(z_val(maxZ_coord-1),1),:);
        if(any(dataInThisDepth,"all"))
            r = dataInThisDepth(1,13)*0.04;  
            pos = [dataInThisDepth(1,1),dataInThisDepth(1,2)]; 
            t=0:0.001:(2*pi);  
            t=[t,0];
            BoundaryPlot_handle = plot(pos(1)+r*sin(t),pos(2)+r*cos(t),'k');

            legend([BoundaryPlot_handle(1)],"Eddy Boundary");
        end
        hold off
        
        xlim([minX-2, maxX+2]);
        ylim([minY-2, maxY+2]);

        figure,

        imagesc(x_val, y_val, ow_val(:,:,maxZ_coord)');
        cb = colorbar();
        caxis(colorRange);
        hold on

        quiver(x_val, y_val, u_val(:,:,maxZ_coord)', v_val(:,:,maxZ_coord)', 5,'r');
        daspect([1,1,1])
        hold on
        title("bottom layer velocity field in frame " + num2str(frameIndex));



        dataInThisDepth = data(round(data(:,5),1) == round(z_val(maxZ_coord),1),:);
        if(any(dataInThisDepth,"all"))
            r = dataInThisDepth(1,13)*0.04;  
            pos = [dataInThisDepth(1,1),dataInThisDepth(1,2)]; 
            t=0:0.001:(2*pi);  
            t=[t,0];
            BoundaryPlot_handle = plot(pos(1)+r*sin(t),pos(2)+r*cos(t),'k');

            legend([BoundaryPlot_handle(1)],"Eddy Boundary");
        end
        hold off
        xlabel("longitude");
        ylabel("latitude");

        xlim([minX-2, maxX+2]);
        ylim([minY-2, maxY+2]);
        
        figure,

        imagesc(x_val, y_val, ow_val(:,:,maxZ_coord+1)');
        cb = colorbar();
        caxis(colorRange);
        hold on

        quiver(x_val, y_val, u_val(:,:,maxZ_coord+1)', v_val(:,:,maxZ_coord+1)', 5,'r');
        daspect([1,1,1])
        
        hold on
        
        dataInThisDepth = data(round(data(:,5),1) == round(z_val(maxZ_coord+1),1),:);

        if(any(dataInThisDepth,"all"))
            r = dataInThisDepth(1,13)*0.04;  
            pos = [dataInThisDepth(1,1),dataInThisDepth(1,2)]; 
            t=0:0.001:(2*pi);  
            t=[t,0];
            BoundaryPlot_handle = plot(pos(1)+r*sin(t),pos(2)+r*cos(t),'k');

            legend([BoundaryPlot_handle(1)],"Eddy Boundary");
        end
        hold off
        xlabel("longitude");
        ylabel("latitude");
        title("bottom + 1 layer lvelocity field in frame " + num2str(frameIndex));
        

        
        xlim([minX-2, maxX+2]);
        ylim([minY-2, maxY+2]);

        figure,

        imagesc(x_val, y_val, ow_val(:,:,maxZ_coord+2)');
        cb = colorbar();
        caxis(colorRange);
        hold on
        quiver(x_val, y_val, u_val(:,:,maxZ_coord+2)', v_val(:,:,maxZ_coord+2)', 5,'r');
        daspect([1,1,1])
        
        hold on
        
        
        dataInThisDepth = data(round(data(:,5),1) == round(z_val(maxZ_coord+2),1),:);

        if(any(dataInThisDepth, "all"))
            r = dataInThisDepth(1,13)*0.04;  
            pos = [dataInThisDepth(1,1),dataInThisDepth(1,2)]; 
            t=0:0.001:(2*pi);  
            t=[t,0];
            BoundaryPlot_handle = plot(pos(1)+r*sin(t),pos(2)+r*cos(t),'k');

            legend([BoundaryPlot_handle(1)],"Eddy Boundary");
        end
        hold off
        xlabel("longitude");
        ylabel("latitude");
        
        title("bottom + 2 layer velocity field in frame " + num2str(frameIndex));



        
        xlim([minX-2, maxX+2]);
        ylim([minY-2, maxY+2]);



    end

end

