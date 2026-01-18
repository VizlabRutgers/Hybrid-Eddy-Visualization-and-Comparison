function [] = eddyBioVis_surface(srcData, property, dataFilePath, allEddy, eddyHistory,timeIndex)
%EDDYSUDDENCHANGEVIS 此处显示有关此函数的摘要
%   此处显示详细说明

close all;

x_val = double(ncread(srcData, property.x));
y_val = double(ncread(srcData, property.y));
z_val = double(ncread(srcData, property.z));

suddenChangeThresh = 0;

% figure,
% heatmap(averEccentricityOverTimesteps');
% title("Average Eccentricity of eddies on each layer (depth) over time steps for Red Sea dataset");
% xlabel("Timesteps");
% ylabel("Depth");

currTimestepEddy = allEddy(allEddy(:,15)==timeIndex,:);
close all;
figure,
ax_isosurface_prev = axes();
daspect(ax_isosurface_prev,[1 1 1]);

%%
temp = double(ncread(srcData, property.temp, [1,1,1,1], [length(x_val),length(y_val),1,1]));

boundaries = bwboundaries(isnan(temp));

for k = 1:length(boundaries)
    boundary = boundaries{k};
    y = x_val(boundary(:,1));  
    x = y_val(boundary(:,2));  

    % 使用 fill 进行可视化
    fill(ax_isosurface_prev,y, x, [0.55,0.45,0.25], 'EdgeColor', 'k', 'FaceAlpha', 0.5);
    hold on
end


%%

u_val_prev = ncread(srcData, property.u, [1,1,1,timeIndex],[length(x_val),length(y_val),length(z_val),1]);
v_val_prev = ncread(srcData, property.v, [1,1,1,timeIndex],[length(x_val),length(y_val),length(z_val),1]);
temp_val_prev = ncread(srcData, property.temp, [1,1,1,timeIndex],[length(x_val),length(y_val),length(z_val),1]);
ow_val_prev = calcOW(length(x_val), length(y_val), u_val_prev, v_val_prev);


for caseIndex = 1:1:size(currTimestepEddy,1)
% for caseIndex = 1:1:size(suddenChangeTimesteps,1)
    
    currentEddy = currTimestepEddy(caseIndex,:);

    prevTimestep = currentEddy(15);
    prevObjIndex = currentEddy(16);
    rotationFlag = currentEddy(14);

    prevEddyData = readEddyData(dataFilePath, rotationFlag, prevTimestep,prevObjIndex);

    
    if(round(prevEddyData(1,5),4)~=round(z_val(1),4))
        continue
    end
    
%     prevEddyData = load(dataFilePath+"Seperated Structures/Frame_"+num2str(prevTimestep)+"_eddy_"+num2str(prevEddyIndex)+"_statistic.uocd");
%     nextEddyData = load(dataFilePath+"Seperated Structures/Frame_"+num2str(nextTimestep)+"_eddy_"+num2str(nextEddyIndex)+"_statistic.uocd");
    
    prevEddy_depth = max(prevEddyData(:,5));

    
    
    
    
    [~,minX_index_prev] = ismember(round(min(prevEddyData(:,3)),4), round(x_val,4));
    [~,maxX_index_prev] = ismember(round(max(prevEddyData(:,3)),4), round(x_val,4));
    [~,minY_index_prev] = ismember(round(min(prevEddyData(:,4)),4), round(y_val,4));
    [~,maxY_index_prev] = ismember(round(max(prevEddyData(:,4)),4), round(y_val,4));

    
    

    %% seabed, land and sea surface
    
%     stretchMode = 2;
%     eddyBioVis_inner(prevEddyData,ax_isosurface_prev,z_val,stretchMode,...
%         1,srcData, property);
%     hold on
% 
% 
% 
%     set(ax_isosurface_prev,'ZDir','reverse');
%     set(ax_isosurface_prev, 'YDir', 'normal');
%     view(3)
%     daspect(ax_isosurface_prev,[1,1,100])    
% 
%     zlim(ax_isosurface_prev,[0,300]);
%     cmap = jet(256);       % 原始 colormap（你也可以用 jet、hot 等）
%     grayness = 0.7;           % 0=原始颜色, 1=完全灰色
%     
%     % 混合到灰色（灰色RGB值是[0.5 0.5 0.5]）
%     grayish_cmap = (cmap + 1) /2;    
%     colormap(ax_isosurface_prev,grayish_cmap);
% 

%     grid on
%     title(ax_isosurface_prev,"winding angle eddy detection results over 10 timesteps");
%    
% 
%     lightangle(ax_isosurface_prev,45,-45);
%     lighting(ax_isosurface_prev, 'gouraud');
%     hold off
%     
%     set(ax_isosurface_prev, "FontSize", 22);


    
    
    %% surface radius ring belts
    
%     u_val_next = ncread(srcData, property.u, [1,1,1,nextTimestep],[length(x_val),length(y_val),length(z_val),1]);
%     v_val_next = ncread(srcData, property.v, [1,1,1,nextTimestep],[length(x_val),length(y_val),length(z_val),1]);
%     temp_val_next = ncread(srcData, property.temp, [1,1,1,nextTimestep],[length(x_val),length(y_val),length(z_val),1]);
%     ow_val_next = calcOW(length(x_val), length(y_val), u_val_next, v_val_next);
    
    prevEddyData_surface = prevEddyData(round(prevEddyData(:,5),4) == round(z_val(1),4),:);

    [~,x_surfaceCenter] = ismember(round(prevEddyData_surface(1,1),4), round(x_val,4));
    [~,y_surfaceCenter] = ismember(round(prevEddyData_surface(1,2),4), round(y_val,4));

    [~,x_pts] = ismember(round(prevEddyData_surface(:,3),4), round(x_val,4));
    radius = (max(x_pts) - min(x_pts))/2;
    radius_real = (max(prevEddyData_surface(:,3)) - min(prevEddyData_surface(:,3)))/2;
    radius_list = [ceil(radius/3), radius-ceil(radius/3), ceil(radius)];
    radius_real_list = [radius_real/3, radius_real-radius_real/3, radius_real];


    %circle
    
    [x_grid,y_grid] = meshgrid(1:length(x_val), 1:length(y_val));

    mask_circle = ((x_grid-x_surfaceCenter).^2 + (y_grid-y_surfaceCenter).^2) <= radius_list(1)^2;  
    temp_val_surface = temp_val_prev(:,:,1)';
    temp_inside_mean = mean(temp_val_surface(mask_circle)); 
    temp_result = NaN(size(temp_val_surface));
    temp_result(mask_circle) = temp_inside_mean;
    
    mask_middleRing = (((x_grid-x_surfaceCenter).^2 + (y_grid-y_surfaceCenter).^2) > radius_list(1)^2) & (((x_grid-x_surfaceCenter).^2 + (y_grid-y_surfaceCenter).^2) <= radius_list(2)^2);  
    temp_val_surface = temp_val_prev(:,:,1)';
    temp_middle_mean = mean(temp_val_surface(mask_middleRing)); 
    temp_result(mask_middleRing) = temp_middle_mean;


    mask_outerRing = (((x_grid-x_surfaceCenter).^2 + (y_grid-y_surfaceCenter).^2) > radius_list(2)^2) & (((x_grid-x_surfaceCenter).^2 + (y_grid-y_surfaceCenter).^2) <= radius_list(3)^2);  
    temp_val_surface = temp_val_prev(:,:,1)';
    temp_outer_mean = mean(temp_val_surface(mask_outerRing)); 
    temp_result(mask_outerRing) = temp_outer_mean;

    vals = [0, temp_inside_mean-temp_middle_mean, temp_inside_mean-temp_outer_mean];  % 每个环的值
    radii = [0 radius_real_list(1); radius_real_list(1) radius_real_list(2); radius_real_list(2) radius_real_list(3)]; % 每一行是 [inner, outer]

    hold on;
    theta = linspace(0, 2*pi, 200);
    for i = 1:size(radii,1)
        r1 = radii(i,1);
        r2 = radii(i,2);

        x_out = r2 * cos(theta)+prevEddyData_surface(1,1);
        y_out = r2 * sin(theta)+prevEddyData_surface(1,2);
        z_out = z_val(1) * ones(size(theta));
        x_in = fliplr(r1 * cos(theta))+prevEddyData_surface(1,1);
        y_in = fliplr(r1 * sin(theta))+prevEddyData_surface(1,2);
        z_in = z_val(1) * ones(size(theta));
       

        x_all = [x_out, x_in];
        y_all = [y_out, y_in];
        z_all = [z_out, z_in];

        fill3(ax_isosurface_prev,x_all, y_all, z_all,vals(i), 'EdgeColor', 'none');
        hold on;
    end
    

%     set(ringFig, 'AlphaData', ~isnan(temp_result));
%     axis on;
%     caxis(ax_isosurface_prev,[min(temp_inside_mean, temp_outer_mean),max(temp_inside_mean, temp_outer_mean)]);
%     colormap("parula");



    %% Isosurface




    prev_x = prevEddyData(:,3);
    prev_y = prevEddyData(:,4);
    prev_z = prevEddyData(:,5);
    prev_temp = prevEddyData(:,11);


    bound = boundary(prev_x,prev_y,prev_z,0.9);
    edgeAlpha=0;
    localOb = trisurf(bound,prev_x,prev_y,prev_z,prev_temp, 'FaceAlpha', '0.3', 'EdgeAlpha', 0,'Parent',ax_isosurface_prev);

    [gridY, gridX, gridZ] = meshgrid(y_val(minY_index_prev-10:maxY_index_prev+10), x_val(minX_index_prev-10:maxX_index_prev+10), z_val);

    % prev timestep
    cmap = lines(10);



    voxel_threshold = 100;
    value_threshold = 26.5;
    color = cmap(1,:);
    alpha = 0.5;
    owIndividualIsosurface(ax_isosurface_prev, gridX, gridY, gridZ, minX_index_prev, maxX_index_prev, minY_index_prev, maxY_index_prev, temp_val_prev, voxel_threshold,min(temp_inside_mean, temp_outer_mean), color,alpha);

%     voxel_threshold = 100;
%     value_threshold = -0.05;
%     color = cmap(2,:);
%     alpha = 0.15;
%     owIndividualIsosurface(ax, gridX, gridY, gridZ, minX_index, maxX_index, minY_index, maxY_index, ow_val_prev, voxel_threshold,value_threshold, color,alpha);

    value_threshold = 26.8;
    color = cmap(3,:);
    alpha = 0.9;
    owIndividualIsosurface(ax_isosurface_prev, gridX, gridY, gridZ, minX_index_prev, maxX_index_prev, minY_index_prev, maxY_index_prev, temp_val_prev, voxel_threshold,max(temp_inside_mean, temp_outer_mean), color,alpha);

    prevEddyCenter = sortrows(unique(prevEddyData(:,[1,2,5]),"rows"),3);
    scatter3(prevEddyCenter(:,2),prevEddyCenter(:,1),prevEddyCenter(:,3),50,".k","DisplayName","eddy centers");
    hold on
    plot3(prevEddyCenter(:,2),prevEddyCenter(:,1),prevEddyCenter(:,3),'color',[1,0,0,0.5],'LineWidth',3,"DisplayName","eddy centerline");
    daspect([1 1 300]);    
    hold off;
    title(ax_isosurface_prev, "Isosurface of OW in timestep " + num2str(prevTimestep));
    axis on;
    legend show;

% %     next timestep
%     
%     value_threshold = -0.01;
%     color = cmap(1,:);
%     alpha = 0.1;
%     owIndividualIsosurface(ax_isosurface_next, gridX, gridY, gridZ, minX_index, maxX_index, minY_index, maxY_index, ow_val_next, voxel_threshold,value_threshold, color,alpha);
% 
% %     voxel_threshold = 100;
% %     value_threshold = -0.05;
% %     color = cmap(2,:);
% %     alpha = 0.15;
% %     owIndividualIsosurface(ax, gridX, gridY, gridZ, minX_index, maxX_index, minY_index, maxY_index, ow_val_prev, voxel_threshold,value_threshold, color,alpha);
% 
%     value_threshold = -0.15;
%     color = cmap(3,:);
%     alpha = 0.4;
%     owIndividualIsosurface(ax_isosurface_next, gridX, gridY, gridZ, minX_index, maxX_index, minY_index, maxY_index, ow_val_next, voxel_threshold,value_threshold, color,alpha);
%     
%     nextEddyCenter = sortrows(unique(nextEddyData(:,[1,2,5]),"rows"),3);
%     scatter3(nextEddyCenter(:,2),nextEddyCenter(:,1),nextEddyCenter(:,3),50,".k","DisplayName","eddy centers");
%     hold on
%     plot3(nextEddyCenter(:,2),nextEddyCenter(:,1),nextEddyCenter(:,3),'color',[1,0,0,0.5],'LineWidth',3,"DisplayName","eddy centerline");
%     daspect([1 1 300]);    
%     hold off;
%     title(ax_isosurface_prev, "Isosurface of OW in timestep " + num2str(prevTimestep));
%     legend show;
    
%     hold off;
%     title(ax_isosurface_next, "Isosurface of OW in timestep " + num2str(nextTimestep));
%     legend show;

    set(ax_isosurface_prev, "FontSize", 22);
%     set(ax_isosurface_next, "FontSize", 22);





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
%     ax_prevVertical_vel=subplot(1,4,1);
%     uimagesc(ax_prevVertical_vel, y_val, z_val, squeeze(sqrt(u_val_prev(prevEddy_centerX,:,:).^2+v_val_prev(prevEddy_centerX,:,:).^2))');
%     xlim(ax_prevVertical_vel,[min(prevEddyData(:,2)) - 1.5, max(prevEddyData(:,2)) + 1.5]);
% %     ylim([0,prevEddy_depth]);
%     xlabel("latitude");
%     ylabel("depth (meter)")
%     cb = colorbar();
%     cb.Label.String = "velocity magnitude";
%     title("Vertical profile of velocity magnitude in timestep " + num2str(prevTimestep));
%     daspect([1 300 1]);
%     set(ax_prevVertical_vel, "FontSize", 14);
% 
%     ax_nextVertical_vel=subplot(1,4,2);
%     uimagesc(ax_nextVertical_vel, y_val, z_val, squeeze(sqrt(u_val_next(nextEddy_centerX,:,:).^2+v_val_next(nextEddy_centerX,:,:).^2))');
%     xlim(ax_nextVertical_vel, [min(nextEddyData(:,2)) - 1.5, max(nextEddyData(:,2)) + 1.5]);
% %     ylim([0,prevEddy_depth]);
%     xlabel("latitude");
%     ylabel("depth (meter)")
%     cb = colorbar();
%     cb.Label.String = "velocity magnitude";
%     title("Vertical profile of velocity magnitude in timestep " + num2str(nextTimestep));
%     daspect([1 300 1]);
%     set(ax_nextVertical_vel, "FontSize", 14);
% 
% 
%     ax_prevVertical_OW=subplot(1,4,3);
%     uimagesc(ax_prevVertical_OW, y_val, z_val, squeeze(ow_val_prev(prevEddy_centerX,:,:))');
%     xlim(ax_prevVertical_OW,[min(prevEddyData(:,2)) - 1.5, max(prevEddyData(:,2)) + 1.5]);
% %     ylim([0,prevEddy_depth]);
%     xlabel("latitude");
%     ylabel("depth (meter)")
%     cb = colorbar();
%     cb.Label.String = "OW";
%     caxis([-0.03,0]);
%     title("Vertical profile of OW in timestep " + num2str(prevTimestep));
%     daspect([1 300 1]);
%     set(ax_prevVertical_OW, "FontSize", 14);
% 
%     ax_nextVertical_OW=subplot(1,4,4);
%     uimagesc(ax_nextVertical_OW, y_val, z_val, squeeze(ow_val_next(nextEddy_centerX,:,:))');
%     xlim(ax_nextVertical_OW, [min(nextEddyData(:,2)) - 1.5, max(nextEddyData(:,2)) + 1.5]);
% %     ylim([0,prevEddy_depth]);
%     xlabel("latitude");
%     ylabel("depth (meter)")
%     cb = colorbar();
%     cb.Label.String = "OW";
%     caxis([-0.03,0]);
%     title("Vertical profile of OW in timestep " + num2str(nextTimestep));
%     daspect([1 300 1]);
%     set(ax_nextVertical_OW, "FontSize", 14);
% 
%     %% Horizontal Vis
%     currentLayer = cutoffDepthIndex;
%     suddenChangeHorizontalVis(srcData, property, cutoffDepthIndex, ...
%         prevEddyData,nextEddyData, currentLayer, prevTimestep, nextTimestep, cutoffDepthFlag, ow_val_prev,ow_val_next);
%     hold off
%     
% 
%     currentLayer = cutoffDepthIndex+1;
%     suddenChangeHorizontalVis(srcData, property, cutoffDepthIndex, ...
%         prevEddyData,nextEddyData, currentLayer, prevTimestep, nextTimestep, cutoffDepthFlag, ow_val_prev,ow_val_next);
%     hold off
end
cb = colorbar(ax_isosurface_prev);
cb.Label.String = "Temperature Difference";
caxis(ax_isosurface_prev,[25.5,27]);
xlabel(ax_isosurface_prev,"Longitude");
ylabel(ax_isosurface_prev,"Latitude");
zlabel(ax_isosurface_prev,"Depth");
caxis([-0.5,0.5]);
daspect([1 1 1])
xlim([x_val(1), x_val(end)]);
ylim([y_val(1), y_val(end)]);
ax_isosurface_prev.FontSize=24;
end

