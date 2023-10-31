function [outputArg1,outputArg2] = eddyMethodComparison(owDataFilePath,windingAngleDataFilePath,srcData,eddyPathHistory,frameIndex)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



srcpath=srcData;

x_val = double(ncread(srcpath, "XC"));
y_val = double(ncread(srcpath, "YC"));
z_val = double(ncread(srcpath, "Z_MIT40"));
totalFrames = double(ncread(srcpath, 'T_AX'));




%% ----------Background---------
f1 = figure();
f1.WindowState = 'maximized';
ax1 = axes(f1);

% imagesc(ax1, x_val, y_val, (ow_val(:,:,1)'));
% colorbar();
% clim([-0.02,0]);



v = VideoWriter('eddyMethodComparison_2D_redsea.mp4','MPEG-4');
v.FrameRate=1;
v.Quality=100;
open(v);

%% ----------Loop of timeframes---------
% You could loop here instead of one single frame
% for i = 1:1:length(totalFrames)
for frame =frameIndex


    startLoc = [1,1,1,frame];
    count = [length(x_val),length(y_val),length(z_val),1];
    stride = [1,1,1,1];
    u_val = ncread(srcData, "U", startLoc, count, stride);
    v_val = ncread(srcData, "V",startLoc, count, stride);

    % w_val = ncread(srcpath, "W",startLoc, count, stride);
    % eta_val = ncread(srcpath, "ETA",[1,1,frame], [length(x_val),length(y_val),1], [1,1,1]);
    % temp_val = ncread(srcpath, "TEMP",startLoc, count, stride);
    % salinity_val = ncread(srcpath, "SALT",startLoc, count, stride);
    vel_mag = sqrt(u_val.^2+v_val.^2);
    
    xSize = length(x_val);
    ySize = length(y_val);
    
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
    % Normalization
    ow_val = ow_val./(vel_mag.^2);
    
    [x_Rgrid2D,y_Rgrid2D] = ndgrid(x_val,y_val);
    DataSize_2D = [length(x_val), length(y_val)];
    
    cla(ax1);

    I0 = imagesc(ax1, x_val, y_val, ~(vel_mag(:,:,1)'));
    color_map=[1 1 1; 0.3 0.3 0.3];
    colormap(color_map);
    hold on

    quiver(ax1,x_Rgrid2D,y_Rgrid2D,u_val(:,:,1),v_val(:,:,1),"AutoScaleFactor",1.5,"LineWidth",1);
    set(ax1,'YDir','normal');
    hold on
    
    
    
    
    %% ----------Winding Angle Method---------
    
    
    load(windingAngleDataFilePath, 'eddy_center_x');
    load(windingAngleDataFilePath, 'eddy_center_y');
    load(windingAngleDataFilePath, 'eddy_ellipse_theta');
    load(windingAngleDataFilePath, 'eddy_length_x');
    load(windingAngleDataFilePath, 'eddy_length_y');
    
%     eddy_sizeFilterThresh = 0.06;
%     eddy_longAxis = max(eddy_length_x,eddy_length_x);
%     eddy_filteredFlag = find(eddy_longAxis<eddy_sizeFilterThresh);
%     
%     eddy_center_x(eddy_filteredFlag) = [];
%     eddy_center_y(eddy_filteredFlag) = [];
%     eddy_ellipse_theta(eddy_filteredFlag) = [];
%     eddy_length_x(eddy_filteredFlag) = [];
%     eddy_length_y(eddy_filteredFlag) = [];
%     eddy_longAxis(eddy_filteredFlag) = [];
    
    WindingAngle_Mask = false(DataSize_2D);
    
    for i = 1:1:length(eddy_center_x)
        l1 = eddy_length_x(i);
        l2 = eddy_length_y(i);
        ellipse_theta = eddy_ellipse_theta(i);
        it = 1:361;
        ellipse_x = l1*cosd(it);
        ellipse_y = l2*sind(it);
        % Create rotation matrix
        R = [cosd(ellipse_theta) -sind(ellipse_theta); sind(ellipse_theta) cosd(ellipse_theta)];
        % Rotate points
        eddypts = R*[ellipse_x ; ellipse_y];
        ellipse_x = eddypts(1,:);
        ellipse_y = eddypts(2,:);
        
    %     m_plot(eddy_center_x(i),eddy_center_y(i),...
    %         'kx','linewidth',4,'markersize',12)
    %     hold on
        windingAnlgeMethod = plot(ax1,eddy_center_x(i) + ellipse_x,...
            eddy_center_y(i) + ellipse_y,'-k','LineWidth',4);
        hold on
        % convert coordinates from pixel coordinate system to degree
        % coordinate system
        windingAnglepolyMask(:,:,i) = poly2mask(coordConvert(eddy_center_x(i) + ellipse_x, x_val,3,0.04),coordConvert(eddy_center_y(i) + ellipse_y, y_val, 3,0.04),DataSize_2D(1),DataSize_2D(2));
     
    end
    
    for i = 1:1:length(eddy_center_x)
        WindingAngle_Mask = WindingAngle_Mask | windingAnglepolyMask(:,:,i);
    end
    
    % figure,imshow(WindingAngle_Mask);
    
    hold on
    
    %% ----------Hybrid Method---------
    data = eddyPathHistory(eddyPathHistory(:,15)==frame,:);
    centerX = data(:,1);
    centerY = data(:,2);
    radius = data(:,13);
    [dataRows,~] = size(data);


    for k = 1:dataRows
        thisEddy = data(k,:);
        r = thisEddy(13)*0.04;  
        pos = [thisEddy(1),thisEddy(2)]; 
        t=0:0.001:(2*pi);  
        t=[t,0];
        BoundaryPlot_handle = plot(ax1,pos(1)+r*sin(t),pos(2)+r*cos(t),'r','linewidth',4);
        hold on
    end
    
    hold on 
    %% ----------OW Method---------
    data = load(owDataFilePath + "red_sea_"+num2str(frame)+"_Clean.uocd");
    data(data(:,4)~=z_val(1),:)=[];
    x = data(:,2);
    y = data(:,3);
    
    OWMask = false(DataSize_2D);
    OWMask(sub2ind(DataSize_2D, coordConvert(y, y_val,2,0.04), coordConvert(x, x_val,2)))=true;
%     OW_WindingAngleIOU = sum(OWMask&WindingAngle_Mask)/sum(OWMask|WindingAngle_Mask);
    % figure,imshow(OWMask);
    
    OWBoundary = bwboundaries(OWMask);
    
    for k = 1:length(OWBoundary)
        boundary = OWBoundary{k};
        plot(x_val(boundary(:,2)), y_val(boundary(:,1)), "Color", [0.9290 0.6940 0.1250], 'LineWidth', 4);
        plot(x_val(boundary(:,2)), y_val(boundary(:,1)), "Color", 'g', 'LineWidth', 4);
        hold on
    end
    xlim(ax1,[43,50]);
    ylim(ax1,[10,15]);

    %% Create Video
%     daspect([1 1 500]);
%     F1=getframe(f1);
%     writeVideo(v,F1);
end

end