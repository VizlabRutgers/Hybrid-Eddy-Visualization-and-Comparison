function [] = eddyGlobalVis(allPath, FrameNum, srcData,sizeLimit, property,dataFilePath)
% allPath: All of eddy's path (used for containg every eddy here)
% FrameNum: Decide which frame to be counted in the global statistic
% srcData: The src data of .nc file
% sizeLimit: Filter out the eddy lower than this limit
% patchSize: Decide the size of patch
%   Detailed explanation goes here
%------------------------------------
%     close all;

    % Read data
    x_val = double(ncread(srcData, property.x));
    y_val = double(ncread(srcData, property.y));
    z_val = double(ncread(srcData, property.z));
    startLoc = [1,1,1,FrameNum];
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
%     salinity_val = ncread(srcData,"salt",startLoc, count, stride);

    
    % Get the coast boundary
%     coastRegion=~isnan(temp_val(:,:,1));
%     [B,L] = bwboundaries(coastRegion,'holes');
    
%     velocity_mag_val = (sqrt(u_val.^2 +v_val.^2));
%     
%     
%     u_val_test = u_val(1:xSize-1,:,:) - u_val(2:xSize,:,:);
%     u_x = u_val;
%     u_x(1:xSize-1,:,:) = u_val_test;
%     u_val_test = u_val(:,1:ySize-1,:) - u_val(:,2:ySize,:);
%     u_y = u_val;
%     u_y(:,1:ySize-1,:) = u_val_test;
%     
%     v_val_test = v_val(1:xSize-1,:,:) - v_val(2:xSize,:,:);
%     v_x = v_val;
%     v_x(1:xSize-1,:,:) = v_val_test;
%     v_val_test = v_val(:,1:ySize-1,:) - v_val(:,2:ySize,:);
%     v_y = v_val;
%     v_y(:,1:ySize-1,:) = v_val_test;
%     
%     s_n = (u_x-v_y).^2;
%     s_s = (u_y+v_x).^2;
%     omega = (v_x-u_y).^2;
%     
%     s_n(s_n==0) = NaN;
%     s_s(s_n==0) = NaN;
%     omega(s_n==0) = NaN;
%     ow_val = s_s+s_n-omega;
%     ow_val = ow_val./(velocity_mag_val.^2);
    


    
    % Tranverse every frame
    for frameIndex=FrameNum
        eddyPath=allPath(allPath(:,15)==frameIndex & allPath(:,13)>=sizeLimit.lower,:);
        
%         Code below is used to see statistic of each eddy without patch
%         setting
%         ----------------------------------------------------------
        figure,
        ax1 = axes;


%         uimage(ax1,x_val,y_val,temp_val(:,:,1)',"CDataMapping","scaled");
        hold on
        [x_Rgrid2D,y_Rgrid2D] = ndgrid(x_val(startLoc(1):1:(startLoc(1)+count(1))-1),y_val(startLoc(2):1:(startLoc(2)+count(2))-1));
        % 
        quiver(ax1,x_Rgrid2D,y_Rgrid2D,u_val(:,:,1),v_val(:,:,1),10,'r','AutoScaleFactor',0.5);


        title("Eddy Size for frame"+num2str(1));
        ylabel('Longitude');
        xlabel('Latitude');
        daspect([1,1,1]);
        hold on
        scatter(ax1,eddyPath(:,1),eddyPath(:,2),[],eddyPath(:,13),"filled");
        cb=colorbar();
        cb.Label.String='Eddy Size';
        caxis([4,10]);
        set(ax1,'YDir','normal');

    
        figure,ax2 = axes;
        daspect(ax2, [1 1 500]);
        
        figure,ax3 = axes;
        daspect(ax3, [1 1 500]);
        
        for eddyIndex = 1:1:size(eddyPath,1)


            if(eddyPath(eddyIndex,14)==1)
                data = load(dataFilePath+"Seperated Structures/clockwise/Frame_"+num2str(FrameNum)+"_eddy_"+num2str(eddyPath(eddyIndex,16))+"_statistic.uocd");
            elseif(eddyPath(eddyIndex,14)==0)
                data = load(dataFilePath+"Seperated Structures/counterclockwise/Frame_"+num2str(FrameNum)+"_eddy_"+num2str(eddyPath(eddyIndex,16))+"_statistic.uocd");
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
            
            bound = boundary(x,y,z,0.9);
            localOb = trisurf(bound,x,y,z,sqrt(u_obj.^2+v_obj.^2),'EdgeColor','none', 'FaceAlpha', '0.9','Parent', ax2);
%             localOb.SpecularExponent = 200;
%             localOb.AmbientStrength = 0.8;
            hold(ax2,'on');
            
            scatter(ax3, x,y,z);
            hold(ax3, 'on');

            
        end



%         for i = 1:1:length(size(eddyPath,1))
%             thisEddy=eddyPath(i,:);
%             r = thisEddy(13)*0.04;  
%             pos = [thisEddy(1),thisEddy(2)]; 
%             t=0:0.001:(2*pi);  
%             t=[t,0];
%             BoundaryPlot_handle = plot(pos(1)+r*sin(t),pos(2)+r*cos(t));
%             hold on
%             legend([BoundaryPlot_handle(1)],"Eddy Boundary");
% 
% 
%         end

%         patch(pos(1)+r*sin(t),pos(2)+r*cos(t), h*ones(size(t)),'k');
%         ----------------------------------------------------------

    
    end
                
    set(ax2, 'ZDir', 'reverse');
    set(ax2, 'YDir', 'normal');
    set(ax2, 'XDir', 'reverse');
    daspect(ax2, [1,1,500]);
    hold on
%     quiver(ax2, x_Rgrid2D,y_Rgrid2D,u_val(:,:,1),v_val(:,:,1),10,'r','AutoScaleFactor',0.5);


%     cb = colorbar();
%     cb.Label.String='Average Eddy Size Times Nums';
% %     set(h3, "alphadata",~isnan(eddyOWDivision)');
%     pbaspect([1,1,1]);
%     title("Avereage Eddy Size Times Nums for whole dataset("+int2str(patchSize*7)+"km x "+int2str(patchSize*7)+"km per patch)");
%     ylabel('Lattude');
%     xlabel('Longitude');
% %     caxis([0,0.006]);
%     set(ax7,'YDir','normal');
end

