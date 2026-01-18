function [outputArg1,outputArg2] = ellipseInterpolationVis(srcData, property, eddyIndex, eddyPathHistory,dataFilePath, interpolationRecord, stride4Test, stretchMode)
%UNTITLED 


interpolationRecord_filtered = interpolationRecord([interpolationRecord.eddyPathNumber] == eddyIndex);
interpolationRecord_filtered = interpolationRecord_filtered([interpolationRecord_filtered.interpolationStride] == stride4Test);
% interpolationRecord_filtered = interpolationRecord_filtered([interpolationRecord_filtered.startIndex] == 1);



% interpolationStartFrame = min(cell2mat(interpolationRecord(:,3)));
% interpolationEndFrame = max(cell2mat(interpolationRecord(:,5)));
% 
% interpolationRecord = interpolationRecord(cell2mat(interpolationRecord(:,2)) == stride4Test, :);
% 
% interpolationFrameIndex = interpolationStartFrame:stride4Test:interpolationEndFrame;
% interpolationRecord_filtered = interpolationRecord(ismember(cell2mat(interpolationRecord(:,3)),interpolationFrameIndex),:);


x_val = double(ncread(srcData, property.x));
y_val = double(ncread(srcData, property.y));
z_val = double(ncread(srcData, property.z));





[x_Rgrid3D,y_Rgrid3D,z_Rgrid3D] = ndgrid(x_val,y_val,z_val);


fh1 = figure();
ax1=axes(fh1);
fh2 = figure();
ax2=axes(fh2);

ax1.FontSize=24;
ax2.FontSize=24;

camlight(ax2);
lighting(ax2, 'flat');


hold(ax2,'on');
set(ax2,'ZDir','reverse');
view(3)




%% groundtruth

for differentEddyIndex=eddyIndex
    

    % Read eddy history
    eddyHistory = eddyPathHistory{differentEddyIndex};
    eddyHistory=sortrows(eddyHistory,15);

    
    % Extract the sub-graph of the eddy path

    
    
    
    minX=inf;
    maxX=0;
    minY=inf;
    maxY=0;
    traceHistory=[];
    
    
    
    minX=min(eddyHistory(:,1));
    maxX=max(eddyHistory(:,1));
    minY=min(eddyHistory(:,2));
    maxY=max(eddyHistory(:,2));
    
    
    % for index=1:1:size(eddyHistory,1)
    for index=1:1:10
        % Get each eddy in this eddy path history at this frame 
        % Plot the boundary of this eddy in 2D


        data = load(dataFilePath+"Seperated Structures/Frame_"+num2str(eddyHistory(index,15))+"_eddy_"+num2str(eddyHistory(index,16))+"_statistic.uocd");


        eddyIndividual3DVis_inner(data,ax2,z_val, stretchMode, index, srcData, property);
        hold on
    end      
end




set(ax2,'ZDir','reverse');
set(ax2, 'YDir', 'normal');
view(3)
daspect([1,1,300])    
ylim([-2,2.5]);
xlim([-4,26]);
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

xticks(0:6:54);
xticklabels({'1','2','3','4','5','6','7','8','9','10'});

% grid on;
% set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
% cb.Label.String = "Velocity Magnitude";


% Remove y axis numbering




% lightangle(0,30);
lightangle(45,-45);
lighting(ax2, 'gouraud');
% lighting(ax2, 'flat');

% delete(findall(gcf, 'Type', 'light'));

%% Interpolation

startFrame = interpolationRecord_filtered.startIndex;


for interpolatedEddyIndex = 1:1:size(interpolationRecord_filtered,1)
% for interpolatedEddyIndex = 1
    interpolation_CenterX = interpolationRecord_filtered.interp_data.centerX;
    interpolation_CenterY = interpolationRecord_filtered.interp_data.centerY;
    interpolation_theta =interpolationRecord_filtered.interp_data.theta;
    interpolation_LengthX = interpolationRecord_filtered.interp_data.lengthX;
    interpolation_LengthY = interpolationRecord_filtered.interp_data.lengthY;


    % Center interpolation

    if(stretchMode == 1)
        interpolatedSurfaceOffset_X = interpolation_CenterX(1)*9;
        interpolatedSurfaceOffset_Y = interpolation_CenterY(1)*9;    
    else
        interpolatedSurfaceOffset_X = (cell2mat(interpolationRecord_filtered(interpolatedEddyIndex,4)) - startFrame)*3;
        interpolatedSurfaceOffset_Y = (cell2mat(interpolationRecord_filtered(interpolatedEddyIndex,4)) - startFrame)*3;   
    end



%     scatter(ax2,interpolationCenter_X(1) + interpolatedSurfaceOffset_X, interpolationCenter_Y(1) + interpolatedSurfaceOffset_Y, 25,'k',"filled");
    interpolation_center_handle = plot3(ax2,interpolation_CenterX(:) + interpolatedSurfaceOffset_X, interpolation_CenterY(:) + interpolatedSurfaceOffset_Y, z_val(1:1:length(interpolation_CenterX)),'k','LineWidth',3);
    hold on

    % Radius interpolation

    for depthIndex = 1:1:length(interpolation_CenterX)

        eddy_ellipse_theta = interpolation_theta(depthIndex);
        eddy_length_x = interpolation_LengthX(depthIndex);
        eddy_length_y = interpolation_LengthY(depthIndex);
        eddy_center_x = interpolation_CenterX(depthIndex) + interpolatedSurfaceOffset_X;
        eddy_center_y = interpolation_CenterY(depthIndex) + interpolatedSurfaceOffset_Y;
        

        %     d(:,1);
        %     d(:,2);
        it = 1:360;
        ellipse_x_shift = eddy_length_x*cosd(it);
        ellipse_y_shift = eddy_length_y*sind(it);
        % Create rotation matrix
        R = [cosd(eddy_ellipse_theta) -sind(eddy_ellipse_theta); sind(eddy_ellipse_theta) cosd(eddy_ellipse_theta)];
        % Rotate points
        eddypts = R*[ellipse_x_shift ; ellipse_y_shift];
        ellipse_x_shift = eddypts(1,:);
        ellipse_y_shift = eddypts(2,:);
        eddy_elipse_x = eddy_center_x + ellipse_x_shift;
        eddy_elipse_y = eddy_center_y + ellipse_y_shift;
        
        plot3(ax2,eddy_elipse_x,eddy_elipse_y,z_val(depthIndex*ones(360,1)),'r.')
        hold on
    end



end

xticklabels("")
yticklabels("")
zticklabels("")
% Set major Grid lines
ax.GridLineStyle = '-';
ax.GridColor = 'k';
ax.GridAlpha = 1;
grid on;
% Set minor Grid lines
ax.MinorGridLineStyle = '-';
ax.MinorGridColor = 'b';
ax.MinorGridAlpha = 0.5;
grid minor;

legend([interpolation_BoundaryPlot_handle,interpolation_center_handle], "interpolation result", "interpolation center");
end

