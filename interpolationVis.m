function [outputArg1,outputArg2] = interpolationVis(srcData, property, eddyIndex, eddyPathHistory,dataFilePath, interpolationRecord, stride4Test, stretchMode)
%UNTITLED 



interpolationStartFrame = min(cell2mat(interpolationRecord(:,3)));
interpolationEndFrame = max(cell2mat(interpolationRecord(:,5)));

interpolationRecord = interpolationRecord(cell2mat(interpolationRecord(:,2)) == stride4Test, :);

interpolationFrameIndex = interpolationStartFrame:stride4Test:interpolationEndFrame;
interpolationRecord_filtered = interpolationRecord(ismember(cell2mat(interpolationRecord(:,3)),interpolationFrameIndex),:);


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
    
    
    for index=1:1:size(eddyHistory,1)
        % Get each eddy in this eddy path history at this frame 
        % Plot the boundary of this eddy in 2D

        if(eddyHistory(index,14)==1)
            data = load(dataFilePath+"Seperated Structures/clockwise/Frame_"+num2str(eddyHistory(index,15))+"_eddy_"+num2str(eddyHistory(index,16))+"_statistic.uocd");
        elseif(eddyHistory(index,14)==0)
            data = load(dataFilePath+"Seperated Structures/counterclockwise/Frame_"+num2str(eddyHistory(index,15))+"_eddy_"+num2str(eddyHistory(index,16))+"_statistic.uocd");
        else
            error('Error: Can not find corresponding eddy data');
        end

        thisEddy = eddyHistory(index,:);
        r = thisEddy(13)*abs(x_val(2) - x_val(1));  
        pos = [thisEddy(1),thisEddy(2)]; 
        t=0:0.001:(2*pi);  
        t=[t,0];
        BoundaryPlot_handle = plot(ax1,pos(1)+r*sin(t),pos(2)+r*cos(t),'k','linewidth',6);
        hold(ax1,'on');

        eddyIndividual3DVis_inner(data,ax2,z_val, stretchMode, index)
    end      
end

set(ax2,'ZDir','reverse');
view(3)
daspect([1 1 500])

cb = colorbar();


% grid on;
% set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
% cb.Label.String = "Velocity Magnitude";


% Remove y axis numbering


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

camlight(ax2);
lighting(ax2, 'flat');

%% Interpolation

startFrame = min(cell2mat(interpolationRecord_filtered(:,3)));


for interpolatedEddyIndex = 1:1:size(interpolationRecord_filtered,1)
    interpolationCenter_X = cell2mat(interpolationRecord_filtered(interpolatedEddyIndex,18));
    interpolationCenter_Y = cell2mat(interpolationRecord_filtered(interpolatedEddyIndex,19));


    % Center interpolation

    if(stretchMode == 1)
        interpolatedSurfaceOffset_X = interpolationCenter_X(1)*9;
        interpolatedSurfaceOffset_Y = interpolationCenter_Y(1)*9;    
    else
        interpolatedSurfaceOffset_X = (cell2mat(interpolationRecord_filtered(interpolatedEddyIndex,4)) - startFrame)*3;
        interpolatedSurfaceOffset_Y = (cell2mat(interpolationRecord_filtered(interpolatedEddyIndex,4)) - startFrame)*3;   
    end



    scatter(ax2,interpolationCenter_X(1) + interpolatedSurfaceOffset_X, interpolationCenter_Y(1) + interpolatedSurfaceOffset_Y, 25,'k',"filled");
    interpolation_center_handle = plot3(ax2,interpolationCenter_X(:) + interpolatedSurfaceOffset_X, interpolationCenter_Y(:) + interpolatedSurfaceOffset_Y, z_val(1:1:length(interpolationCenter_X)),'k','LineWidth',3);


    % Radius interpolation

    for depthIndex = 1:1:length(interpolationCenter_X)

        interpolationRadius = cell2mat(interpolationRecord_filtered(interpolatedEddyIndex,15));
        
        pos = [interpolationCenter_X(depthIndex) + interpolatedSurfaceOffset_X,interpolationCenter_Y(depthIndex) + interpolatedSurfaceOffset_Y]; 
        t=0:0.001:(2*pi);  
        t=[t,0];
        depth = ones(1,length(t))*z_val(depthIndex);
        interpolation_BoundaryPlot_handle = plot3(ax2,pos(1)+interpolationRadius(depthIndex)*sin(t),pos(2)+interpolationRadius(depthIndex)*cos(t), depth,'r','linewidth',3);
        hold(ax2,'on');
    end



end

legend([interpolation_BoundaryPlot_handle,interpolation_center_handle], "interpolation result", "interpolation center");
end

