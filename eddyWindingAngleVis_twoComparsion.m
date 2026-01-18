function [] = eddyWindingAngleVis_twoComparsion(srcData, property, dataFilePath, dataFilePath_new,allEddy, allEddy_new, eddyHistory, eddyHistory_new, eddyPathIndex, eddyPathIndex_new)
%EDDYSUDDENCHANGEVIS 此处显示有关此函数的摘要
%   此处显示详细说明

% close all;

x_val = double(ncread(srcData, property.x));
y_val = double(ncread(srcData, property.y));
z_val = double(ncread(srcData, property.z));
time_val = ncread(srcData, property.time);

% figure,
% heatmap(averEccentricityOverTimesteps');
% title("Average Eccentricity of eddies on each layer (depth) over time steps for Red Sea dataset");
% xlabel("Timesteps");
% ylabel("Depth");





for timestep = 1:1:length(time_val)
%     close all;   
    %% Isosurface
    u_val_old = ncread(srcData, property.u, [1,1,1,timestep],[length(x_val),length(y_val),length(z_val),1]);
    v_val_old = ncread(srcData, property.v, [1,1,1,timestep],[length(x_val),length(y_val),length(z_val),1]);
    ow_val_old = calcOW(length(x_val), length(y_val), u_val_old, v_val_old);
    [gridX, gridY, gridZ] = meshgrid(y_val, x_val, z_val);

%     figure
%     quiver(x_val, y_val, u_val_old, v_val_old);
%     
    
    
    % old timestep
    cmap = lines(10);
    figure,
    ax_isosurface_old = axes();
    voxel_threshold = 100;
    value_threshold = -0.01;
    color = cmap(1,:);
    alpha = 0.1;
    
    minX_index = 1;
    maxX_index = length(x_val);
    minY_index = 1;
    maxY_index = length(y_val);
    
    
    ow_high_handle = owIndividualIsosurface(ax_isosurface_old, minX_index, maxX_index, minY_index, maxY_index, ow_val_old, voxel_threshold,value_threshold, color,alpha, x_val, y_val, z_val);

%     voxel_threshold = 100;
%     value_threshold = -0.05;
%     color = cmap(2,:);
%     alpha = 0.15;
%     owIndividualIsosurface(ax, gridX, gridY, gridZ, minX_index, maxX_index, minY_index, maxY_index, ow_val_old, voxel_threshold,value_threshold, color,alpha);

    value_threshold = -0.2;
    color = cmap(3,:);
    alpha = 0.3;
    ow_low_handle = owIndividualIsosurface(ax_isosurface_old, minX_index, maxX_index, minY_index, maxY_index, ow_val_old, voxel_threshold,value_threshold, color,alpha, x_val, y_val, z_val);
    
    xlim([43.5, 45.4]);
    ylim([10.5 13]);
    zlim([1,1000]);
    view(3);
    
    %% Individual eddies
    for oldEddyIndex = eddyPathIndex
        currEddyPath = eddyHistory{oldEddyIndex};

        oldEddy = currEddyPath(currEddyPath(:,15)==timestep,:);

        if(isempty(oldEddy))
            continue;
        end
        oldEddyObjIndex = oldEddy(16);

        oldEddyData = load(dataFilePath+"Seperated Structures/Frame_"+num2str(timestep)+"_eddy_"+num2str(oldEddyObjIndex)+"_statistic.uocd");

        oldEddy_depth = max(oldEddyData(:,5));

        oldEddyCenter = sortrows(unique(oldEddyData(:,[1,2,5]),"rows"),3);
        scatter3(ax_isosurface_old,oldEddyCenter(:,1),oldEddyCenter(:,2),oldEddyCenter(:,3),50,".k");
        hold on
        oldCenterline_handle = plot3(ax_isosurface_old,oldEddyCenter(:,1),oldEddyCenter(:,2),oldEddyCenter(:,3),'color',[0.5625, 0.0000, 0.5625],'LineWidth',3);  
        hold on;

    end
    
    for newEddyndex = eddyPathIndex_new
        currEddyPath_new = eddyHistory_new{newEddyndex};

        newEddy = currEddyPath_new(currEddyPath_new(:,15)==timestep,:);

        if(isempty(newEddy))
            continue;
        end
        newEddyObjIndex = newEddy(16);

        newEddyData = load(dataFilePath_new+"Seperated Structures/Frame_"+num2str(timestep)+"_eddy_"+num2str(newEddyObjIndex)+"_statistic.uocd");

        newEddy_depth = max(newEddyData(:,5));

        newEddyCenter = sortrows(unique(newEddyData(:,[1,2,5]),"rows"),3);
        scatter3(ax_isosurface_old,newEddyCenter(:,1),newEddyCenter(:,2),newEddyCenter(:,3),50,".","MarkerFaceColor","#77AC30");
        hold on
        newCenterline_handle = plot3(ax_isosurface_old,newEddyCenter(:,1),newEddyCenter(:,2),newEddyCenter(:,3),'color',"#77AC30",'LineWidth',3);
  
        hold on;
    end
%     cutoffDepth = min(oldEddy_depth,newEddy_depth);
%     if(oldEddy_depth < newEddy_depth)
%         cutoffDepthFlag=-1;
%     elseif(oldEddy_depth > newEddy_depth)
%         cutoffDepthFlag = 1;
%     else
%         cutoffDepthFlag = 0;
%     end
%     
%     [~,cutoffDepthIndex] = ismember(cutoffDepth, z_val);
    
    
%     [~,minX_index_old] = ismember(round(min(oldEddyData(:,3)),4), round(x_val,4));
%     [~,maxX_index_old] = ismember(round(max(oldEddyData(:,3)),4), round(x_val,4));
%     [~,minY_index_old] = ismember(round(min(oldEddyData(:,4)),4), round(y_val,4));
%     [~,maxY_index_old] = ismember(round(max(oldEddyData(:,4)),4), round(y_val,4));
% 
%     [~,minX_index_new] = ismember(round(min(newEddyData(:,3)),4), round(x_val,4));
%     [~,maxX_index_new] = ismember(round(max(newEddyData(:,3)),4), round(x_val,4));
%     [~,minY_index_new] = ismember(round(min(newEddyData(:,4)),4), round(y_val,4));
%     [~,maxY_index_new] = ismember(round(max(newEddyData(:,4)),4), round(y_val,4));
% 
%     minX_index = min(minX_index_old,minX_index_new);
%     maxX_index = max(maxX_index_old,maxX_index_new);
%     minY_index = min(minY_index_old,minY_index_new);
%     maxY_index = max(maxY_index_old,maxY_index_new);




    %% 3D Vis
%     fh2 = figure();
%     ax2=axes(fh2);
%     
%     stretchMode = 2;
%     eddyIndividual3DVis_inner(oldEddyData,ax2,z_val,stretchMode,...
%         1,srcData, property);
%     hold on
%     eddyIndividual3DVis_inner(newEddyData,ax2,z_val,stretchMode,...
%         2,srcData, property)
%     hold on
% 
% 
%     set(ax2,'ZDir','reverse');
%     set(ax2, 'YDir', 'normal');
%     view(3)
%     daspect([1,1,300])    
%     ylim([-2,2.5]);
%     xlim([-4,10]);
%     zlim([0,1000]);
%     colormap(jet(5));
%     cb = colorbar();
%     cb.Label.String = "Velocity Magnitude";
%     caxis([0,0.5])
%     xlabel("Timestep");
%     ylabel("Actual size in Degree");
%     zlabel("Depth");
%     grid on
%     title("winding angle eddy detection results over 10 timesteps");
%     
%     xticks(0:6:6);%     currentLayer = cutoffDepthIndex;
%     suddenChangeHorizontalVis(srcData, property, cutoffDepthIndex, ...
%         oldEddyData,newEddyData, currentLayer, oldTimestep, newTimestep, cutoffDepthFlag, ow_val_old,ow_val_new);
%     hold off
%     
% 
%     currentLayer = cutoffDepthIndex+1;
%     suddenChangeHorizontalVis(srcData, property, cutoffDepthIndex, ...
%         oldEddyData,newEddyData, currentLayer, oldTimestep, newTimestep, cutoffDepthFlag, ow_val_old,ow_val_new);
%     hold off
%     xticklabels({num2str(timestep),num2str(timestep)});
% 
%     lightangle(45,-45);
%     lighting(ax2, 'gouraud');
%     hold off
%     

% 
%     set(ax2, "FontSize", 22);







    
    % new timestep
%     ax_isosurface_new = subplot(1,2,2);
%     value_threshold = -0.01;
%     color = cmap(1,:);
%     alpha = 0.1;
%     owIndividualIsosurface(ax_isosurface_new, gridX, gridY, gridZ, minX_index, maxX_index, minY_index, maxY_index, ow_val_old, voxel_threshold,value_threshold, color,alpha);

%     voxel_threshold = 100;
%     value_threshold = -0.05;
%     color = cmap(2,:);
%     alpha = 0.15;
%     owIndividualIsosurface(ax, gridX, gridY, gridZ, minX_index, maxX_index, minY_index, maxY_index, ow_val_old, voxel_threshold,value_threshold, color,alpha);
% 
%     value_threshold = -0.15;
%     color = cmap(3,:);
%     alpha = 0.4;
%     owIndividualIsosurface(ax_isosurface_new, gridX, gridY, gridZ, minX_index, maxX_index, minY_index, maxY_index, ow_val_old, voxel_threshold,value_threshold, color,alpha);
    
%     hold off;
%     title(ax_isosurface_new, "Isosurface of OW in timestep " + num2str(timestep));
%     legend show;


%     set(ax_isosurface_new, "FontSize", 22);

    %% Centerline 
%     figure,



    %% Volume rendering on OW
%     ow_val_excludeNaN = -ow_val_filtered;
%     mask = ~isfinite(ow_val_excludeNaN);
%     ow_val_excludeNaN(mask) = 0;
% 
%     CC = bwconncomp(ow_val_excludeNaN, 26);
%     % 获取每个区域的体素数
%     stats = regionprops3(CC, 'Volume');  % 体素个数
%     
%     % 设置体积阈值（小于该值的区域将被删除）
%     min_voxel_count = 500;
%     
%     % 找出大的区域
%     large_idx = find(stats.Volume >= min_voxel_count);
%     
%     % 创建新 mask，只保留大的区域
%     cleaned_mask = ismember(labelmatrix(CC), large_idx);
%     
% 
% 
%     % 应用 mask 到原始数据
%     data_cleaned = ow_val_excludeNaN;
%     data_cleaned(~cleaned_mask) = 0;
%     h = volumeViewer(data_cleaned);


    

    %% Vertical profile
%     figure,
%     ax_oldVertical_vel=subplot(1,4,1);
%     uimagesc(ax_oldVertical_vel, y_val, z_val, squeeze(sqrt(u_val_old(oldEddy_centerX,:,:).^2+v_val_old(oldEddy_centerX,:,:).^2))');
%     xlim(ax_oldVertical_vel,[min(oldEddyData(:,2)) - 1.5, max(oldEddyData(:,2)) + 1.5]);
% %     ylim([0,oldEddy_depth]);
%     xlabel("latitude");
%     ylabel("depth (meter)")
%     cb = colorbar();
%     cb.Label.String = "velocity magnitude";
%     title("Vertical profile of velocity magnitude in timestep " + num2str(oldTimestep));
%     daspect([1 300 1]);
%     set(ax_oldVertical_vel, "FontSize", 14);
% 
%     ax_newVertical_vel=subplot(1,4,2);
%     uimagesc(ax_newVertical_vel, y_val, z_val, squeeze(sqrt(u_val_old(newEddy_centerX,:,:).^2+v_val_old(newEddy_centerX,:,:).^2))');
%     xlim(ax_newVertical_vel, [min(newEddyData(:,2)) - 1.5, max(newEddyData(:,2)) + 1.5]);
% %     ylim([0,oldEddy_depth]);
%     xlabel("latitude");
%     ylabel("depth (meter)")
%     cb = colorbar();
%     cb.Label.String = "velocity magnitude";
%     title("Vertical profile of velocity magnitude in timestep " + num2str(newTimestep));
%     daspect([1 300 1]);
%     set(ax_newVertical_vel, "FontSize", 14);
% 
% 
%     ax_oldVertical_OW=subplot(1,4,3);
%     uimagesc(ax_oldVertical_OW, y_val, z_val, squeeze(ow_val_old(oldEddy_centerX,:,:))');
%     xlim(ax_oldVertical_OW,[min(oldEddyData(:,2)) - 1.5, max(oldEddyData(:,2)) + 1.5]);
% %     ylim([0,oldEddy_depth]);
%     xlabel("latitude");
%     ylabel("depth (meter)")
%     cb = colorbar();
%     cb.Label.String = "OW";
%     caxis([-0.03,0]);
%     title("Vertical profile of OW in timestep " + num2str(oldTimestep));
%     daspect([1 300 1]);
%     set(ax_oldVertical_OW, "FontSize", 14);
% 
%     ax_newVertical_OW=subplot(1,4,4);
%     uimagesc(ax_newVertical_OW, y_val, z_val, squeeze(ow_val_old(newEddy_centerX,:,:))');
%     xlim(ax_newVertical_OW, [min(newEddyData(:,2)) - 1.5, max(newEddyData(:,2)) + 1.5]);
% %     ylim([0,oldEddy_depth]);
%     xlabel("latitude");
%     ylabel("depth (meter)")
%     cb = colorbar();
%     cb.Label.String = "OW";
%     caxis([-0.03,0]);
%     title("Vertical profile of OW in timestep " + num2str(newTimestep));
%     daspect([1 300 1]);
%     set(ax_newVertical_OW, "FontSize", 14);

    %% Horizontal Vis
%     currentLayer = cutoffDepthIndex;
%     suddenChangeHorizontalVis(srcData, property, cutoffDepthIndex, ...
%         oldEddyData,newEddyData, currentLayer, oldTimestep, newTimestep, cutoffDepthFlag, ow_val_old,ow_val_new);
%     hold off
%     
% 
%     currentLayer = cutoffDepthIndex+1;
%     suddenChangeHorizontalVis(srcData, property, cutoffDepthIndex, ...
%         oldEddyData,newEddyData, currentLayer, oldTimestep, newTimestep, cutoffDepthFlag, ow_val_old,ow_val_new);
%     hold off


    %% others
%     title(ax_isosurface_old, "Isosurface of OW in timestep " + num2str(timestep));
%     legend show;
%     title(ax_isosurface_old, "Isosurface of OW in timestep " + num2str(timestep));
    set(ax_isosurface_old, "FontSize", 22);
    daspect([1 1 300]);  
    set(gca, "YDir", "normal");
%     legend([ow_low_handle, ow_high_handle,oldCenterline_handle], ...
%         "OW isosurface with value at" + num2str(-0.2), "OW isosurface with value at" + num2str(-0.01),...
%         "original winding angle centerline");
    legend([ow_low_handle, ow_high_handle,oldCenterline_handle,newCenterline_handle], ...
        "OW isosurface with value at" + num2str(-0.2), "OW isosurface with value at" + num2str(-0.01),...
        "original winding angle centerline","new winding angle centerline");
%     legend([ow_low_handle, ow_high_handle,oldCenterline_handle], ...
%         "OW isosurface with value at" + num2str(-0.2), "OW isosurface with value at" + num2str(-0.01),...
%         "winding angle eddy centerline");
end



end

