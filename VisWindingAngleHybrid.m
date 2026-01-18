function [] = VisWindingAngleHybrid(srcData, property, dataFilePath, allEddy, windingAngleSummary,windingAngleTracks, matches, unMatchHybrid, timeIndex,eddy_ellipse_objects)
%VISWINDINGANGLEHYBRID 此处显示有关此函数的摘要
%   此处显示详细说明

close all;
eddySummary = windingAngleSummary;


x_val = double(ncread(srcData, property.x));
y_val = double(ncread(srcData, property.y));
z_val = double(ncread(srcData, property.z));




fh1 = figure();
fh1.WindowState = 'maximized';
ax1=axes(fh1);
set(gca,'FontSize',16);

coastFaces = load('coastFaceData_RS.mat','coastFaces').coastFaces;
coastVerts = load('coastVertData_RS.mat','coastVerts').coastVerts;
seabed = patch('Faces',coastFaces, 'Vertices',coastVerts, ...  
'FaceColor',[.8, .8, .8],...
'EdgeColor','none',...
'FaceAlpha',0.5, ...
'Parent', ax1);
hold on

eddySummaryFrame = eddySummary{timeIndex};
for depthIndex = 1:1:length(eddySummaryFrame)
    eddySummaryOnLayer = eddySummaryFrame{depthIndex};
    if(isempty(eddySummaryOnLayer))
        continue;
    else
        scatter3(ax1,eddySummaryOnLayer{1}, eddySummaryOnLayer{2},z_val(depthIndex)*ones(1,size(eddySummaryOnLayer{1},2)),20, "filled");
        hold on
    end
end


for trackIndex = 1:length(eddy_ellipse_objects)
    currentTrack = windingAngleTracks{trackIndex};
    currentTrackCoord = cell2mat(currentTrack(1:3,:));
    plot3(ax1,currentTrackCoord(1,:), currentTrackCoord(2,:), z_val(currentTrackCoord(3,:)), 'k','linewidth',2);
    
    eddy_ellipse_pts = eddy_ellipse_objects{trackIndex};
    scatter3(eddy_ellipse_pts(:,1),eddy_ellipse_pts(:,2),eddy_ellipse_pts(:,3), "b.");
    
    hold on
end



hybridEddyInFrame = allEddy(allEddy(:,15)==timeIndex,:);
hybridEddyInFrameMatched = hybridEddyInFrame(ismember(hybridEddyInFrame(:,16)+1, matches(:,2)),:);



for index=1:1:size(hybridEddyInFrameMatched,1)
    % Get each eddy in this eddy path history at this frame 
    % Plot the boundary of this eddy in 2D

    if(hybridEddyInFrameMatched(index,14)==1)
        data = load(dataFilePath+"Seperated Structures/clockwise/Frame_"+num2str(hybridEddyInFrameMatched(index,15))+"_eddy_"+num2str(hybridEddyInFrameMatched(index,16))+"_statistic.uocd");
    elseif(hybridEddyInFrameMatched(index,14)==0)
        data = load(dataFilePath+"Seperated Structures/counterclockwise/Frame_"+num2str(hybridEddyInFrameMatched(index,15))+"_eddy_"+num2str(hybridEddyInFrameMatched(index,16))+"_statistic.uocd");
    else
        error('Error: Can not find corresponding eddy data');
    end

    for depthIndex = 1:1:length(unique(data(:,5)))

        eddyOnLayer_gd = data(data(:,5) == z_val(depthIndex),:);
        centerX = mean(eddyOnLayer_gd(:,3));
        centerY = mean(eddyOnLayer_gd(:,4));

        pos = [centerX,centerY]; 
        t=0:0.001:(2*pi);  
        t=[t,0];
        radius =(max(eddyOnLayer_gd(:,3)) - min(eddyOnLayer_gd(:,3)))/2;  
        depth = ones(1,length(t))*z_val(depthIndex);
        hybridBoundaryPlot_handle = plot3(ax1,pos(1)+radius*sin(t),pos(2)+radius*cos(t), depth,'r','linewidth',1.5);
        hold('on');
    end
end    


if(~isempty(unMatchHybrid))
    hybridEddyInFrameNotMatched = hybridEddyInFrame(ismember(hybridEddyInFrame(:,16)+1, unMatchHybrid),:);
    for index=1:1:size(hybridEddyInFrameNotMatched,1)
        % Get each eddy in this eddy path history at this frame 
        % Plot the boundary of this eddy in 2D

        if(hybridEddyInFrameNotMatched(index,14)==1)
            data = load(dataFilePath+"Seperated Structures/clockwise/Frame_"+num2str(hybridEddyInFrameNotMatched(index,15))+"_eddy_"+num2str(hybridEddyInFrameNotMatched(index,16))+"_statistic.uocd");
        elseif(hybridEddyInFrameNotMatched(index,14)==0)
            data = load(dataFilePath+"Seperated Structures/counterclockwise/Frame_"+num2str(hybridEddyInFrameNotMatched(index,15))+"_eddy_"+num2str(hybridEddyInFrameNotMatched(index,16))+"_statistic.uocd");
        else
            error('Error: Can not find corresponding eddy data');
        end

        for depthIndex = 1:1:length(unique(data(:,5)))

            eddyOnLayer_gd = data(data(:,5) == z_val(depthIndex),:);
            centerX = mean(eddyOnLayer_gd(:,3));
            centerY = mean(eddyOnLayer_gd(:,4));

            pos = [centerX,centerY]; 
            t=0:0.001:(2*pi);  
            t=[t,0];
            radius =(max(eddyOnLayer_gd(:,3)) - min(eddyOnLayer_gd(:,3)))/2;  
            depth = ones(1,length(t))*z_val(depthIndex);
            hybridBoundaryPlot_NotMatched_handle = plot3(ax1,pos(1)+radius*sin(t),pos(2)+radius*cos(t), depth,'b','linewidth',1.5);
            hold('on');
        end
    end    
    legend([hybridBoundaryPlot_handle,hybridBoundaryPlot_NotMatched_handle, seabed], "Matched Hybrid Eddy Boundary", "Matched Hybrid Eddy Boundary", "seabed",'location', 'northeast');
else
    legend([hybridBoundaryPlot_handle, seabed], "Matched Hybrid Eddy Boundary", "seabed", 'location', 'northeast');
end




cb = colorbar();
set(ax1, "ZDir", "reverse");
cb.Label.String = "Eccentricity";
xlabel(ax1,"longitude");
ylabel(ax1,"latitude");
zlabel(ax1,"depth(layer)");
daspect([1 1 500]);
set(gca, 'CameraPosition', [36,-122,-12130]);
grid on;

title("Frame" + num2str(timeIndex));
hold off




end

