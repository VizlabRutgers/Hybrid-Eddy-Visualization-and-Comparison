function [distanceSum, distanceRatio,areaSum,areaRatio, statMatrix] = eddyHorizontalMethodComparison_Stat(owDataFilePath,srcData,...
    eddySummary_oldWindingAngle, eddySummary_newWindingAngle,allEddy_hybrid,...
    frameIndex, property,depth)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% close all

srcpath=srcData;

x_val = double(ncread(srcpath, property.x));
y_val = double(ncread(srcpath, property.y));
z_val = double(ncread(srcpath, property.z));
totalFrames = double(ncread(srcpath, property.time));
[x_grid, y_grid] = meshgrid(x_val,y_val);





%% ----------Background---------
f1 = figure();
f1.WindowState = 'maximized';
ax1 = axes(f1);

% imagesc(ax1, x_val, y_val, (ow_val(:,:,1)'));
% colorbar();
% clim([-0.02,0]);



% v = VideoWriter('eddyMethodComparison_2D_redsea.mp4');
% v.FrameRate=1;
% v.Quality=100;
% open(v);

%% ----------Loop of timeframes---------
% You could loop here instead of one single frame
% for i = 1:1:length(totalFrames)
distanceSum = [];
distanceRatio = [];
areaSum = [];
areaRatio = [];
centerOld_X=[];
centerOld_Y=[];
centerNew_X=[];
centerNew_Y=[];
timestep=[];
effectiveRadius=[];
originalArea=[];
 
matchSum = 0;
for frame = frameIndex
    

    startLoc = [1,1,depth,frame];
    count = [length(x_val),length(y_val),1,1];
    stride = [1,1,1,1];
    u_val = ncread(srcData, property.u, startLoc, count, stride);
    v_val = ncread(srcData, property.v,startLoc, count, stride);

    % w_val = ncread(srcpath, "W",startLoc, count, stride);
    % eta_val = ncread(srcpath, "ETA",[1,1,frame], [length(x_val),length(y_val),1], [1,1,1]);
    % temp_val = ncread(srcpath, "TEMP",startLoc, count, stride);
    % salinity_val = ncread(srcpath, "SALT",startLoc, count, stride);
    vel_mag = sqrt(u_val.^2+v_val.^2);
    
    xSize = length(x_val);
    ySize = length(y_val);
    
    ow_val = calcOW(length(x_val), length(y_val), u_val, v_val);
    


    x_test = x_val(1:length(x_val));
    y_test = y_val(1:length(y_val));

    u_test = u_val(1:length(x_val),1:length(y_val),1);
    v_test = v_val(1:length(x_val),1:length(y_val),1);

    u_test(isnan(u_test)) = 0;
    v_test(isnan(v_test)) = 0;

    vel_mag_test = sqrt(u_test.^2+v_test.^2);



    [x_Rgrid2D,y_Rgrid2D] = ndgrid(x_val,y_val);
    DataSize_2D = [length(x_val), length(y_val)];
    
    cla(ax1);

%     I0 = uimagesc(ax1, x_test, y_test, (~(vel_mag_test')));
    background = pcolor(ax1, x_grid, y_grid, double(vel_mag_test'==0));
    background.FaceColor=[0.3 0.3 0.3];
    background.EdgeColor = 'none';
    alpha(background, double(vel_mag_test'==0)); 
    
    hold on


    % color quiver
    q = quiver(ax1,x_Rgrid2D(1:3:end,1:3:end),y_Rgrid2D(1:3:end,1:3:end),u_test(1:3:end,1:3:end),v_test(1:3:end,1:3:end), "AutoScaleFactor",1.5,"LineWidth",1);
    set(ax1,'YDir','normal');
    hold on
    
    x1 = xlabel(ax1,'Logitude');
    y1 = ylabel(ax1,'Latitude');
    daspect(ax1,[1,1,1]);

    ax1.FontSize=24;
    x1.FontSize = 48;
    y1.FontSize = 48;

    
    
    %% Colored Quiver
    %// Compute the magnitude of the vectors
%     mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), ...
%                 reshape(q.WData, numel(q.UData), [])).^2, 2));
%     
%     cmap = jet(255);
% 
%     %// Get the current colormap
%     currentColormap = cmap;
%     colormap(cmap);
%     
%     %// Now determine the color to make each arrow using a colormap
%     [~, ~, ind] = histcounts(mags, size(currentColormap, 1));
%     
%     %// Now map this to a colormap to get RGB
%     cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
%     cmap(:,:,4) = 255;
%     cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
%     
%     %// We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
%     set(q.Head, ...
%         'ColorBinding', 'interpolated', ...
%         'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'
%     
%     %// We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
%     set(q.Tail, ...
%         'ColorBinding', 'interpolated', ...
%         'ColorData', reshape(cmap(1:2,:,:), [], 4).');
% 
%     cb = colorbar();
% %     colormap(ax2,"hot");
%     cb.Label.String="Velocity Magnitude";
%     c    distanceSum_inFrame = 0;
%     axis(ax1,[0,0.3]);
%     cb.FontSize = 28;
%  
%% ----------OW Background---------
%     [ow_contour, ow_contour_handle] = contour(x_Rgrid2D,y_Rgrid2D, ow_val, [-0.2, -0.2], 'color','#EDB120', "linewidth", 2);
%     contour(x_Rgrid2D,y_Rgrid2D, ow_val', [-0.01, -0.01]);



    %% ----------OW Method---------
%     data = load(owDataFilePath + "red_sea_"+num2str(frame)+"_Clean.uocd");
%     data(data(:,4)~=z_val(depth),:)=[];
%     x = data(:,2);
%     y = data(:,3);
%     
%     OWMask = false(DataSize_2D);
%     OWMask(sub2ind(DataSize_2D, coordConvert(y, y_val,2,0.04), coordConvert(x, x_val,2)))=true;
% %     OW_WindingAngleIOU = sum(OWMask&WindingAngle_Mask)/sum(OWMask|WindingAngle_Mask);
%     % figure,imshow(OWMask);
% %     load("Redsea_eddy_clk_data.mat");
% %     load("Redsea_eddy_conclk_data.mat");
% %     load("Redsea_eddy_data.mat");
% %     load("Redsea_eddy_path.mat");
% %     load("Redsea_eddy_history.mat");
% %     load("Redsea_eddy_graph.mat");
% 
%     OWBoundary = bwboundaries(OWMask);
%     
%     for k = 1:length(OWBoundary)
%         boundary = OWBoundary{k};
%         OW_boundary = plot(x_val(boundary(:,2)), y_val(boundary(:,1)), "Color", [0.9290, 0.6940, 0.1250], 'LineWidth', 4);
%         hold on
%     end
%     xlim(ax1,[43,50]);
%     ylim(ax1,[10,15]);
  
    %% ----------Old Winding Angle Method---------
    eddyHistory_old_currentFrame = eddySummary_oldWindingAngle{frame}{depth};
    currentCenter_X_old = eddyHistory_old_currentFrame{1};
    currentCenter_Y_old = eddyHistory_old_currentFrame{2};
    currentCenter_Z_old = eddyHistory_old_currentFrame{3};
    currentoffset_X_old = eddyHistory_old_currentFrame{4};
    currentoffset_Y_old = eddyHistory_old_currentFrame{5};
    currentCenter_longAxis = eddyHistory_old_currentFrame{6};
    currentCenter_shortAxis = eddyHistory_old_currentFrame{7};
    
    
    for index = 1:1:size(eddyHistory_old_currentFrame{1},2)
        oldWinding_centers = scatter(ax1,currentCenter_X_old(index), currentCenter_Y_old(index),50,"k",'filled');
        hold on
        oldWinding_boundary = plot(ax1, currentCenter_X_old(index)+currentoffset_X_old{index}, currentCenter_Y_old(index)+currentoffset_Y_old{index}, 'color',"k",'linewidth',2);
        area_old(index) = polyarea(currentCenter_X_old(index)+currentoffset_X_old{index}, currentCenter_Y_old(index)+currentoffset_Y_old{index});
        hold on
    end
    
    % figure,imshow(WindingAngle_Mask);
    
    hold on   
%     legend(ax1,[oldWinding_boundary,oldWinding_centers, ow_contour_handle,q], ...
%         "Original winding angle boundary", "Original winding angle center","OW contour with value at" + num2str(-0.2),"Velocity field");

    
    %% ----------New Winding Angle Method---------
    
    eddyHistory_new_currentFrame = eddySummary_newWindingAngle{frame}{depth};
    currentCenter_X_new = eddyHistory_new_currentFrame{1};
    currentCenter_Y_new = eddyHistory_new_currentFrame{2};
    currentCenter_Z_new = eddyHistory_new_currentFrame{3};
    currentboundary_X = eddyHistory_new_currentFrame{6};
    currentboundary_Y = eddyHistory_new_currentFrame{7};
    
    
    for index = 1:1:size(eddyHistory_new_currentFrame{1},2)
        newWinding_centers = scatter(ax1,currentCenter_X_new(index), currentCenter_Y_new(index),50,"r",'filled');
        hold on
        newWinding_boundary = plot(ax1, currentboundary_X{index}, currentboundary_Y{index}, 'color','r','linewidth',2);
        hold on
        area_new(index) = polyarea(currentboundary_X{index}, currentboundary_Y{index});
    end
    
    % figure,imshow(WindingAngle_Mask);
    
    hold on  
    
    %% ----------Hybrid Method---------
%     data = allEddy_hybrid(allEddy_hybrid(:,15)==frame,:);
%     centerX = data(:,1);
%     centerY = data(:,2);
%     radius = data(:,13);
%     [dataRows,~] = size(data);
% 
% 
%     for k = 1:dataRows
%         thisEddy = data(k,:);
%         r = thisEddy(13)*0.04;  
%         pos = [thisEddy(1),thisEddy(2)]; 
%         t=0:0.001:(2*pi);  
%         t=[t,0];
%         hybrid_boundary= plot(ax1,pos(1)+r*sin(t),pos(2)+r*cos(t),'color',[0.4660, 0.6740, 0.1880],'linewidth',4);
%         hold on
%     end
%     
%     hold on 
    
    %% Others
%     legend([OW_boundary, oldWinding_boundary,newWinding_boundary,hybrid_boundary], "OW method", "Old winding angle method", "New winding angle method", "Hybrid method");
    legend([oldWinding_centers,newWinding_centers, oldWinding_boundary,newWinding_boundary,q], ...
        "Original winding angle center", "New winding angle center", ...
        "Original winding angle boundary", "New winding angle boundary", "Velocity field");

    
    %% Statistics
    oldEddies_centerCoord = [currentCenter_X_old',currentCenter_Y_old'];
    newEddies_centerCoord = [currentCenter_X_new',currentCenter_Y_new'];
    
    
    distanceMap = pdist2(newEddies_centerCoord, oldEddies_centerCoord, 'euclidean');
    
    distanceThreshold = 0.15;
    [matches,unMatchCurrent,unMatchPrevious] = matchpairs(distanceMap, distanceThreshold);
    

    for matchIndex = 1:size(matches, 1)

        oldCenter = oldEddies_centerCoord(matches(matchIndex,2),:);
        newCenter = newEddies_centerCoord(matches(matchIndex,1),:);
        
        longAxis = currentCenter_longAxis(matches(matchIndex,1));
        shortAxis = currentCenter_shortAxis(matches(matchIndex,1));

        if(longAxis<0.001 || shortAxis<0.001)
            continue;
        end
        
        distanceSum(end+1) = norm(oldCenter-newCenter);
        distanceRatio(end+1) = norm(oldCenter-newCenter)/sqrt(longAxis*shortAxis);
        centerOld_X(end+1) = currentCenter_X_old(matches(matchIndex,2));
        centerOld_Y(end+1) = currentCenter_Y_old(matches(matchIndex,2));
        centerNew_X(end+1) = currentCenter_X_new(matches(matchIndex,1));
        centerNew_Y(end+1) = currentCenter_Y_new(matches(matchIndex,1));
        timestep(end+1) = frame;
        effectiveRadius(end+1) = sqrt(longAxis*shortAxis);
        originalArea(end+1) = area_old(matches(matchIndex,2));
        
        if(isnan(abs(area_old(matches(matchIndex,2)) - area_new(matches(matchIndex,1)))))
            areaSum(end+1) = -1;
            areaRatio(end+1) = -1;
            matchSum = matchSum + 1;
            continue;
        end
        areaSum(end+1) = abs(area_old(matches(matchIndex,2)) - area_new(matches(matchIndex,1)));
        areaRatio(end+1) = abs(area_old(matches(matchIndex,2)) - area_new(matches(matchIndex,1)))/area_old(matches(matchIndex,2));


        
        matchSum = matchSum + 1;
%         if(matchSum == 15462)
%             matchSum = 15462;
%         end
    end
    
    
    %% Create Video
%     daspect([1 1 500]);
%     F1=getframe(f1);
%     writeVideo(v,F1);
end


fullmap = parula(256);    % 默认 parula
idx = round(linspace(1,256,3));
cmap3 = fullmap(idx, :);  % 取蓝—绿—黄三色

colormap(cmap3);

distanceSum = distanceSum';
distanceRatio = distanceRatio';

areaSum = areaSum';
areaRatio = areaRatio';

statMatrix = [effectiveRadius', distanceRatio, originalArea', areaRatio, ...
    timestep', centerOld_X', centerOld_Y', centerNew_X', centerNew_Y'];

% figure, violinplot(distanceSum);
% title("Boxplot of center distance between two winding angle approaches");
% 
% figure, violinplot(distanceRatio);
% title("Boxplot showing the relative center distance between the two methods (geometry average of original method as the baseline)");
% 
% figure, violinplot(areaSum);
% title("Boxplot of area difference between two winding angle approaches");
% 
% figure, violinplot(areaRatio);
% title("Boxplot showing the relative area difference between the two methods (original method as the baseline)");


end