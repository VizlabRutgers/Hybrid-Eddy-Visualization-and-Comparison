function [] = eddySuddenChangeStat(srcData, property, dataFilePath, allEddy, windingAngleSummary)
%EDDYSUDDENCHANGEVIS 此处显示有关此函数的摘要
%   此处显示详细说明

close all;

eddySummary = windingAngleSummary;

averEccentricityOverTimesteps = zeros(10,50);


x_val = double(ncread(srcData, property.x));
y_val = double(ncread(srcData, property.y));
z_val = double(ncread(srcData, property.z));


for timeIndex=1:1:length(eddySummary)
    eddySummaryFrame = eddySummary{timeIndex};
    averageEccentricityOnLayer = cellfun(@(x) mean(x(end, :)), eddySummaryFrame);
    averEccentricityOverTimesteps(timeIndex,:) = averageEccentricityOnLayer;
end

% figure,
% heatmap(averEccentricityOverTimesteps');
% title("Average Eccentricity of eddies on each layer (depth) over time steps for Red Sea dataset");
% xlabel("Timesteps");
% ylabel("Depth");

for timeIndex=1:1:1
% for timeIndex=1:1:length(eddySummary)
    figure,
    eddySummaryFrame = eddySummary{timeIndex};
    for depthIndex = 1:1:length(eddySummaryFrame)
        eddySummaryOnLayer = eddySummaryFrame{depthIndex};
        if(isempty(eddySummaryOnLayer))
            continue;
        else
            scatter3(eddySummaryOnLayer(1,:), eddySummaryOnLayer(2,:),z_val(depthIndex)*ones(1,length(eddySummaryOnLayer)),20, "filled", "CData",eddySummaryOnLayer(4,:));
            hold on
        end
    end

    hybridEddyInFrame = allEddy(allEddy(:,15)==timeIndex,:);
    
    for index=1:1:size(hybridEddyInFrame,1)
        % Get each eddy in this eddy path history at this frame 
        % Plot the boundary of this eddy in 2D

        if(hybridEddyInFrame(index,14)==1)
            data = load(dataFilePath+"Seperated Structures/clockwise/Frame_"+num2str(hybridEddyInFrame(index,15))+"_eddy_"+num2str(hybridEddyInFrame(index,16))+"_statistic.uocd");
        elseif(hybridEddyInFrame(index,14)==0)
            data = load(dataFilePath+"Seperated Structures/counterclockwise/Frame_"+num2str(hybridEddyInFrame(index,15))+"_eddy_"+num2str(hybridEddyInFrame(index,16))+"_statistic.uocd");
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
            interpolation_BoundaryPlot_handle = plot3(pos(1)+radius*sin(t),pos(2)+radius*cos(t), depth,'r','linewidth',3);
            hold('on');
        end
    end    


    coastFaces = load('coastFaceData_RS.mat','coastFaces').coastFaces;
    coastVerts = load('coastVertData_RS.mat','coastVerts').coastVerts;
    seabed = patch('Faces',coastFaces, 'Vertices',coastVerts, ...  
    'FaceColor',[.8, .8, .8],...
    'EdgeColor','none',...
    'FaceAlpha',0.5);

    cb = colorbar();
    set(gca, "ZDir", "reverse");
    cb.Label.String = "Eccentricity";
    xlabel("longitude");
    ylabel("latitude");
    zlabel("depth(layer)");
    hold off;
    title("Frame" + num2str(timeIndex));
    hold off
end

end

