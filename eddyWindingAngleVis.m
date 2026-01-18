function [] = eddyWindingAngleVis(srcData, property, dataFilePath, allEddy, eddyHistory,eddyPathIndex)
%EDDYSUDDENCHANGEVIS 此处显示有关此函数的摘要
%   此处显示详细说明

close all;

x_val = double(ncread(srcData, property.x));
y_val = double(ncread(srcData, property.y));
z_val = double(ncread(srcData, property.z));

% figure,
% heatmap(averEccentricityOverTimesteps');
% title("Average Eccentricity of eddies on each layer (depth) over time steps for Red Sea dataset");
% xlabel("Timesteps");
% ylabel("Depth");

currEddyPath = eddyHistory{eddyPathIndex};

for i = 1:1:size(currEddyPath,1)
    eddyTimeStepTemp = currEddyPath(i,15);
    eddyIndexTemp = currEddyPath(i,16);
    eddyData = load(dataFilePath+"Seperated Structures/Frame_"+num2str(eddyTimeStepTemp)+"_eddy_"+num2str(eddyIndexTemp)+"_statistic.uocd");
    depthRecord(i,1) = max(eddyData(:,5));
    depthRecord(i,4) = eddyTimeStepTemp;
    [~,depthRecord(i,2)] = ismember(depthRecord(i,1),z_val);
end




for caseIndex = 1:1:size(currEddyPath,1)
    close all;
    currentSuddenChangeCase = currEddyPath(caseIndex,:);
    oldTimestep = currentSuddenChangeCase(15);

    oldEddy = currEddyPath(currEddyPath(:,15)==oldTimestep,:);

    
    [~,oldEddy_centerX] = min(abs(oldEddy(1)-x_val));


    oldEddyIndex = oldEddy(16);

    
    oldEddyData = load(dataFilePath+"Seperated Structures/Frame_"+num2str(oldTimestep)+"_eddy_"+num2str(oldEddyIndex)+"_statistic.uocd");
  
    oldEddy_depth = max(oldEddyData(:,5));
 
    
    
    
    [~,cutoffDepthIndex] = ismember(oldEddy_depth, z_val);
    
    
    [~,minX_index_old] = ismember(round(min(oldEddyData(:,3)),4), round(x_val,4));
    [~,maxX_index_old] = ismember(round(max(oldEddyData(:,3)),4), round(x_val,4));
    [~,minY_index_old] = ismember(round(min(oldEddyData(:,4)),4), round(y_val,4));
    [~,maxY_index_old] = ismember(round(max(oldEddyData(:,4)),4), round(y_val,4));


    minX_index = minX_index_old;
    maxX_index = maxX_index_old;
    minY_index = minY_index_old;
    maxY_index = maxY_index_old;




    %% 3D Vis
    fh2 = figure();
    ax2=axes(fh2);
    
    stretchMode = 2;
    eddyIndividual3DVis_inner(oldEddyData,ax2,z_val,stretchMode,...
        1,srcData, property);
    hold on


    set(ax2,'ZDir','reverse');
    set(ax2, 'YDir', 'normal');
    view(3)
    daspect([1,1,300])    
    ylim([-2,2.5]);
    xlim([-4,10]);
    zlim([0,1000]);
    colormap(jet(5));
    cb = colorbar();
    cb.Label.String = "Velocity Magnitude";
    caxis([0,0.5])
    xlabel("Timestep");
    ylabel("Actual size in Degree");
    zlabel("Depth");
    grid on
    title("winding angle eddy detection results over 10 timesteps");
    
    xticks(0:6:6);
    xticklabels({num2str(oldTimestep)});

    lightangle(45,-45);
    lighting(ax2, 'gouraud');
    hold off
    
    u_val_old = ncread(srcData, property.u, [1,1,1,oldTimestep],[length(x_val),length(y_val),length(z_val),1]);
    v_val_old = ncread(srcData, property.v, [1,1,1,oldTimestep],[length(x_val),length(y_val),length(z_val),1]);
    ow_val_old = calcOW(length(x_val), length(y_val), u_val_old, v_val_old);

    set(ax2, "FontSize", 22);





    %% Isosurface

    [gridX, gridY, gridZ] = meshgrid(y_val(minY_index-10:maxY_index+10), x_val(minX_index-10:maxX_index+10), z_val);

    % old timestep
    cmap = lines(10);
    figure,
    ax_isosurface_old = subplot(1,2,1);
    voxel_threshold = 100;
    value_threshold = -0.01;
    color = cmap(1,:);
    alpha = 0.1;
    owIndividualIsosurface(ax_isosurface_old, gridX, gridY, gridZ, minX_index, maxX_index, minY_index, maxY_index, ow_val_old, voxel_threshold,value_threshold, color,alpha);

%     voxel_threshold = 100;
%     value_threshold = -0.05;
%     color = cmap(2,:);
%     alpha = 0.15;
%     owIndividualIsosurface(ax, gridX, gridY, gridZ, minX_index, maxX_index, minY_index, maxY_index, ow_val_old, voxel_threshold,value_threshold, color,alpha);

    value_threshold = -0.15;
    color = cmap(3,:);
    alpha = 0.3;
    owIndividualIsosurface(ax_isosurface_old, gridX, gridY, gridZ, minX_index, maxX_index, minY_index, maxY_index, ow_val_old, voxel_threshold,value_threshold, color,alpha);

    oldEddyCenter = sortrows(unique(oldEddyData(:,[1,2,5]),"rows"),3);
    scatter3(oldEddyCenter(:,2),oldEddyCenter(:,1),oldEddyCenter(:,3),50,".k","DisplayName","eddy centers");
    hold on
    plot3(oldEddyCenter(:,2),oldEddyCenter(:,1),oldEddyCenter(:,3),'color',[1,0,0,0.5],'LineWidth',3,"DisplayName","eddy centerline");
    daspect([1 1 300]);    
    hold off;
    title(ax_isosurface_old, "Isosurface of OW in timestep " + num2str(oldTimestep));
    legend show;

    
 

    set(ax_isosurface_old, "FontSize", 22);

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
    figure,
    ax_oldVertical_vel=subplot(1,4,1);
    uimagesc(ax_oldVertical_vel, y_val, z_val, squeeze(sqrt(u_val_old(oldEddy_centerX,:,:).^2+v_val_old(oldEddy_centerX,:,:).^2))');
    xlim(ax_oldVertical_vel,[min(oldEddyData(:,2)) - 1.5, max(oldEddyData(:,2)) + 1.5]);
%     ylim([0,oldEddy_depth]);
    xlabel("latitude");
    ylabel("depth (meter)")
    cb = colorbar();
    cb.Label.String = "velocity magnitude";
    title("Vertical profile of velocity magnitude in timestep " + num2str(oldTimestep));
    daspect([1 300 1]);
    set(ax_oldVertical_vel, "FontSize", 14);



    ax_oldVertical_OW=subplot(1,4,3);
    uimagesc(ax_oldVertical_OW, y_val, z_val, squeeze(ow_val_old(oldEddy_centerX,:,:))');
    xlim(ax_oldVertical_OW,[min(oldEddyData(:,2)) - 1.5, max(oldEddyData(:,2)) + 1.5]);
%     ylim([0,oldEddy_depth]);
    xlabel("latitude");
    ylabel("depth (meter)")
    cb = colorbar();
    cb.Label.String = "OW";
    caxis([-0.03,0]);
    title("Vertical profile of OW in timestep " + num2str(oldTimestep));
    daspect([1 300 1]);
    set(ax_oldVertical_OW, "FontSize", 14);

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
end



end

