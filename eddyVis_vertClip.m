function [] = eddyVis_vertClip(pathIndex,G,eddyPath,eddyPathHistory,dataFilePath,srcData, property)
% pathindex: Get Index of individual path history
% eddyPathHistory: The trace history of this eddy
% dataFilePath: The base path datafile for Feature Tracking
% srcData: The src data of .nc file
%   Detailed explanation goes here
%------------------------------------


% Video Recorder

% v2 = VideoWriter('eddyHistory_3D_redsea.mp4');
% v2.FrameRate=1;
% v2.Quality=100;
% open(v2);


% Read coordinates
x_val = double(ncread(srcData, property.x));
y_val = double(ncread(srcData, property.y));
z_val = double(ncread(srcData, property.z));
time_val = double(ncread(srcData, property.time));
resolution = 0.04;

% zVal_Unique = unique(z_val,'rows');
% zVal_even=arrayfun(@(x) find(zVal_Unique==x), z_val);
% zVal_even=zVal_even.*0.1;

CenterHistory=[];
historyCounter=1;

fh1 = figure();
% fh1.WindowState = 'maximized';
ax1=axes(fh1);
% fh2 = figure();
% fh2.WindowState = 'maximized';
% ax2=axes(fh2);

ax1.FontSize=24;
% ax2.FontSize=24;
% 
% camlight(ax2);
% lighting(ax2, 'flat');

pause(1);

% Compute seabed data

% coastRegion = ncread(srcData,property.temp,[1,1,1,1], [length(x_val),length(y_val),length(z_val),1], [1,1,1,1]);
% coast2D=~isnan(coastRegion(:,:,1));
% [B,L] = bwboundaries(coast2D,'holes');
% coastRegion(isnan(coastRegion))=-1;
% [coastFaces,coastVerts]=isosurface(x_val,y_val,z_val,permute(coastRegion,[2,1,3]),-1);
% 
% save('coastFaceData_RS.mat','coastFaces');
% save('coastVertData_RS.mat','coastVerts');

% Load seabed data
% coastFaces = load('coastFaceData_NA.mat','coastFaces').coastFaces;
% coastVerts = load('coastVertData_NA.mat','coastVerts').coastVerts;

% coastFaces = load('coastFaceData_NP.mat','coastFaces').coastFaces;
% coastVerts = load('coastVertData_NP.mat','coastVerts').coastVerts;


maxFrame = cellfun(@(x) max(x(:,15)), eddyPathHistory);
maxFrame = max(maxFrame);

minFrame = cellfun(@(x) min(x(:,15)), eddyPathHistory);
minFrame = min(minFrame);

maxX = cellfun(@(x) max(x(:,1)), eddyPathHistory);
maxX = max(maxX);

minX = cellfun(@(x) min(x(:,1)), eddyPathHistory);
minX = min(minX);

maxY = cellfun(@(x) max(x(:,2)), eddyPathHistory);
maxY = max(maxY);

minY = cellfun(@(x) min(x(:,2)), eddyPathHistory);
minY = min(minY);

radiusRecordInDepth = [];

% Here the i indicates the time frame in eddy history
for i = minFrame:1:2
% for i = 1
    cla(ax1);
%     cla(ax2);

%     if (i == maxFrame)
%         i = maxFrame;
%     end

    FrameNum=i;
    startLoc = [1,1,1,FrameNum];
    count = [length(x_val),length(y_val),1,1];
    stride = [1,1,1,1];
    
    startLoc_2D = [1,1,FrameNum];
    count_2D = [length(x_val),length(y_val),1];
    stride_2D = [1,1,1];


    u_val = ncread(srcData, property.u, startLoc, count, stride);
    v_val = ncread(srcData, property.v,startLoc, count, stride);
%     w_val = ncread(srcData, "W",startLoc, count, stride);
%     eta_val = ncread(srcData,"ETA",startLoc_2D, count_2D, stride_2D);
%     temp_val = ncread(srcData,property.temp,startLoc, count, stride);
%     salinity_val = ncread(srcData,"salt",startLoc, count, stride);
%     [x_Rgrid3D,y_Rgrid3D,z_Rgrid3D] = ndgrid(x_val,y_val,z_val);
%     velocity_mag_val = (sqrt(u_val.^2 +v_val.^2));


    
    % Get quvier map

    



    for differentEddyIndex=1:1:size(eddyPathHistory,1)
        
        

        % Read eddy history
        eddyHistory = eddyPathHistory{differentEddyIndex};
        eddyHistory=sortrows(eddyHistory,15);

        subpath=eddyPath{pathIndex(differentEddyIndex)};
        subGraphMap=subgraph(G,subpath);
        firstEddyNode = int2str(eddyHistory(1,15))+"."+int2str(eddyHistory(1,16)+1);
        firstEddyNodeIndex=find(table2cell(subGraphMap.Nodes)==firstEddyNode);
        

        % Extract the sub-graph of the eddy path

        
        
        
        minX=inf;
        maxX=0;
        minY=inf;
        maxY=0;
        traceHistory=[];
        
    
        % Get the ith frame's eddies (could be more than 1)
        eddyInThisFrame=eddyHistory(eddyHistory(:,15)==i,:);
        eddyInThisFrameNum=size(eddyInThisFrame,1);
        maxRadius=max(eddyHistory(1,13)*0.04);
        
        
        
        

        
        
        for index=1:1:eddyInThisFrameNum
            radiusRecordTemp = [];
            eddyFrameIndex = (eddyHistory(:,15)==i);
            
            if(eddyHistory(eddyFrameIndex,14)==1)
                data = load(dataFilePath+"Seperated Structures/clockwise/Frame_"+num2str(FrameNum)+"_eddy_"+num2str(eddyHistory(eddyFrameIndex,16))+"_statistic.uocd");
            elseif(eddyHistory(eddyFrameIndex,14)==0)
                data = load(dataFilePath+"Seperated Structures/counterclockwise/Frame_"+num2str(FrameNum)+"_eddy_"+num2str(eddyHistory(eddyFrameIndex,16))+"_statistic.uocd");
            else
                error('Error: Can not find corresponding eddy data');
            end

            depth = unique(data(:,5));
            shape_left=[];
            shape_right=[];
            [~, depthIndexTotal] = ismember(round(depth), round(z_val));

            minX=min(data(:,3));
            maxX=max(data(:,3));
            minZ=min(data(:,5));
            maxZ=max(data(:,5));

            if(numel(depthIndexTotal) < 4)
                continue;
            end

            for depthIndex = 1:1:length(depthIndexTotal)
            
                dataInDepth = data(data(:,5)==depth(depthIndex),:);

                x = dataInDepth(:,3);
                y = dataInDepth(:,4);
                z = dataInDepth(:,5);
                ow = dataInDepth(:,6);
                u_obj = dataInDepth(:,7);
                v_obj = dataInDepth(:,8);
                temp = dataInDepth(:,10);
                salt = dataInDepth(:,11);

                radius = (max(x) - min(x))/2;

                shape_left = [shape_left; min(x),depth(depthIndex)];
                shape_right = [shape_right; max(x),depth(depthIndex)];
                radiusRecordTemp = [radiusRecordTemp; radius];
            end


            
            original_data = radiusRecordTemp;                        
            original_index = 1:1:length(radiusRecordTemp);                       
            target_index = 1:1:length(z_val);

            stretched_data = original_data * (numel(target_index) / numel(original_index));
            stretched_index = linspace(min(target_index), max(target_index), numel(original_index));
            radiusRecordTemp = interp1(stretched_index, stretched_data, target_index, 'linear');
            radiusRecordInDepth = [radiusRecordInDepth; radiusRecordTemp];
            
            shape_left(:,1) = (shape_left(:,1) - minX) / (maxX - minX);
            shape_right(:,1) = (shape_right(:,1) - minX) / (maxX - minX);


            shape_left(:,2) = (shape_left(:,2) - minZ) / (maxZ - minZ);
            shape_right(:,2) = (shape_right(:,2) - minZ) / (maxZ - minZ);

%             plot(ax1,shape_left(:,1), shape_left(:,2),"k");
%             hold(ax1,'on');
%             plot(ax1,shape_right(:,1), shape_right(:,2),"k");
%             hold(ax1,'on');
%             plot(ax1, [shape_left(end,1);shape_right(end,1)],[shape_left(end,2);shape_right(end,2)],"k");
%             hold(ax1,'on');
            
            % Get the boundary of this eddy structure
%                 bound = boundary(x,y,z,0.9);
%                 localOb = trisurf(bound,x,y,z,sqrt(u_obj.^2+v_obj.^2),'EdgeColor','none', 'FaceAlpha', '0.9','Parent', ax2);
%                 localOb.SpecularExponent = 200;
%                 localOb.AmbientStrength = 0.8;
%                 hold(ax2,'on');
            

        end
    end



%         for boundaryIndex = 1:1:size(B,1)
%             coastboundary=B{boundaryIndex};
%             plot(ax2,x_val(coastboundary(:,1)), y_val(coastboundary(:,2)),'k','LineWidth',4);
%             hold(ax2,'on');
%         end
end



    % 使用 kmeans 函数进行聚类
    num_classes = 6;
    [idx, centroids] = kmeans(radiusRecordInDepth, num_classes, "Distance","correlation");
    
    for class_index=1:1:num_classes
        [~, lineIndex] = ismember(idx,class_index);
        figure, plot(radiusRecordInDepth(lineIndex==1, :)');
    end



    % 绘制聚类结果
    figure,gscatter(1:size(radiusRecordInDepth,1), 1:size(radiusRecordInDepth,1), idx);
    hold on;
    plot(centroids(:,1), centroids(:,2), 'kx', 'MarkerSize', 15, 'LineWidth', 3);
    legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Centroids');
    xlabel('特征1');
    ylabel('特征2');
    title('K-means 聚类结果');
    hold off;


    
    set(ax1, "YDir", "reverse");
    
%     camlight(ax2);
%     lighting(ax2, 'flat');
%     lightangle(ax2,-37.5,-45);


%     xlim(ax1,[minX-maxRadius-0.5,maxX+maxRadius+0.5]);
%     ylim(ax1,[minY-maxRadius-0.5,maxY+maxRadius+0.5]);

    xlim(ax1,[minX-1,maxX+1]);
    ylim(ax1,[minY-1,maxY+1]);

    if(size(pathIndex)==1)
        title(ax1,"eddy "+int2str(pathIndex)+" from NA dataset Frame"+num2str(FrameNum));
    else
        title(ax1,"eddies from NA dataset Frame"+num2str(FrameNum));
    end
    x1 = xlabel(ax1,'Longitude');
    y1 = ylabel(ax1,'Latitude');
    daspect(ax1,[1,1,1]);

    x1.FontSize = 48;
    y1.FontSize = 48;



    

end



