function [] = eddyVis(pathIndex,G,eddyPath,eddyPathHistory,dataFilePath,srcData, property)
% pathindex: Get Index of individual path history
% eddyPathHistory: The trace history of this eddy
% dataFilePath: The base path datafile for Feature Tracking
% srcData: The src data of .nc file
%   Detailed explanation goes here
%------------------------------------

close all;

% Video Recorder
v = VideoWriter('eddyHistory_2D_redsea.mp4');
v.FrameRate=1;
v.Quality=100;
open(v);

v2 = VideoWriter('eddyHistory_3D_redsea.mp4');
v2.FrameRate=1;
v2.Quality=100;
open(v2);


% Read coordinates
x_val = double(ncread(srcData, property.x));
y_val = double(ncread(srcData, property.y));
z_val = double(ncread(srcData, property.z));
resolution = 0.04;

% zVal_Unique = unique(z_val,'rows');
% zVal_even=arrayfun(@(x) find(zVal_Unique==x), z_val);
% zVal_even=zVal_even.*0.1;

CenterHistory=[];
historyCounter=1;

fh1 = figure();
fh1.WindowState = 'maximized';
ax1=axes(fh1);
fh2 = figure();
fh2.WindowState = 'maximized';
ax2=axes(fh2);

ax1.FontSize=24;
ax2.FontSize=24;

camlight(ax2);
lighting(ax2, 'flat');

pause(5);

% Compute seabed data

% coastRegion = ncread(srcData,property.temp,[1,1,1,1], [length(x_val),length(y_val),length(z_val),1], [1,1,1,1]);
% coastRegion(isnan(coastRegion))=0;
% coast2D=~(coastRegion(:,:,1));
% [B,L] = bwboundaries(coast2D,'holes');
% coastRegion(~logical(coastRegion))=-1;
% [coastFaces,coastVerts]=isosurface(x_val,y_val,z_val,permute(coastRegion,[2,1,3]),-1);
% 
% save('coastFaceData_NA.mat','coastFaces');
% save('coastVertData_NA.mat','coastVerts');

% Load seabed data
coastFaces = load('coastFaceData_NA.mat','coastFaces').coastFaces;
coastVerts = load('coastVertData_NA.mat','coastVerts').coastVerts;

% coastFaces = load('coastFaceData_NP.mat','coastFaces').coastFaces;
% coastVerts = load('coastVertData_NP.mat','coastVerts').coastVerts;

% coastFaces = load('coastFaceData_RS.mat','coastFaces').coastFaces;
% coastVerts = load('coastVertData_RS.mat','coastVerts').coastVerts;

maxFrame = length(ncread(srcData, property.time));


% Here the i indicates the time frame in eddy history
% for i = 1:1:maxFrame
for i = 1:1:26
% for i = 1
    cla(ax1);
    cla(ax2);

%     if (i == maxFrame)
%         i = maxFrame;
%     end

    FrameNum=i;
    startLoc = [1,1,1,FrameNum];
    count = [length(x_val),length(y_val),length(z_val),1];
    stride = [1,1,1,1];
    
    startLoc_2D = [1,1,FrameNum];
    count_2D = [length(x_val),length(y_val),1];
    stride_2D = [1,1,1];


%     u_val = ncread(srcData, property.u, startLoc, count, stride);
%     v_val = ncread(srcData, property.v,startLoc, count, stride);
%     w_val = ncread(srcData, "W",startLoc, count, stride);
%     eta_val = ncread(srcData,"ETA",startLoc_2D, count_2D, stride_2D);
    temp_val = ncread(srcData,property.temp,startLoc, count, stride);
%     salinity_val = ncread(srcData,"salt",startLoc, count, stride);
%     [x_Rgrid3D,y_Rgrid3D,z_Rgrid3D] = ndgrid(x_val,y_val,z_val);
%     velocity_mag_val = (sqrt(u_val.^2 +v_val.^2));
    % Compute OW value
    % xmin = x_val(1);
    % ymin = y_val(1);
    % 
    
    % 
    % 
    % u_val_test = u_val(1:xSize-1,:,:) - u_val(2:xSize,:,:);
    % u_x = u_val;
    % u_x(1:xSize-1,:,:) = u_val_test;
    % u_val_test = u_val(:,1:ySize-1,:) - u_val(:,2:ySize,:);
    % u_y = u_val;
    % u_y(:,1:ySize-1,:) = u_val_test;
    % 
    % v_val_test = v_val(1:xSize-1,:,:) - v_val(2:xSize,:,:);
    % v_x = v_val;
    % v_x(1:xSize-1,:,:) = v_val_test;
    % v_val_test = v_val(:,1:ySize-1,:) - v_val(:,2:ySize,:);
    % v_y = v_val;
    % v_y(:,1:ySize-1,:) = v_val_test;
    % 
    % s_n = (u_x-v_y).^2;
    % s_s = (u_y+v_x).^2;
    % omega = (v_x-u_y).^2;fh2 = figure();

    % 
    % s_n(s_n==0) = NaN;
    % s_s(s_n==0) = NaN;
    % omega(s_n==0) = NaN;
    % ow_val = s_s+s_n-omega;
    % ow_val = ow_val./(velocity_mag_val.^2);

    %% Background (quiver/seabed)
    % Get quvier map
    [x_Rgrid2D,y_Rgrid2D] = ndgrid(x_val(startLoc(1):1:(startLoc(1)+count(1))-1),y_val(startLoc(2):1:(startLoc(2)+count(2))-1));
%     % 
%     quiver(ax1,x_Rgrid2D,y_Rgrid2D,u_val(:,:,1),v_val(:,:,1),'r','AutoScaleFactor',3);

    
%     gray_rgb = repmat(reshape([0.5 0.5 0.5], 1, 1, 3), size(temp_val(:,:,1)'));
% 
%     % gray land
%     imshow(gray_rgb, 'XData', [x_val(1) x_val(end)], 'YData', [y_val(1) y_val(end)]);
%     hold on;

    surface_temp = uimagesc(ax1, x_val,y_val, temp_val(:,:,1)');
    hold(ax1,'on');
    set(ax1, "YDir", "normal");
    daspect(ax1, [1 1 1]);
%     set(surface_temp, 'AlphaData', ~isnan(temp_val(:,:,1)'));
    colormap(jet);
    hold(ax1, 'on');

    % Draw the sea bed
    seabed = patch(ax2,'Faces',coastFaces, 'Vertices',coastVerts, ...  
    'FaceColor',[0.7,0.7,0.7],...
    'EdgeColor','none',...
    'FaceAlpha',0.35, 'DisplayName', 'Seabed');
    isonormals(temp_val, seabed);
    % seabed.SpecularColorReflectance = 0;
    %     seabed.SpecularExponent = 80;
    seabed.AmbientStrength = 0.1;

    daspect(ax2,[1,1,800]);
    
    
    hold(ax2,'on');
    
    
    %% Stream
%     velThresh=0.8;
%     velFlow = (velocity_mag_val-mean(mean(velocity_mag_val,'omitnan'),'omitnan'))./velocity_mag_val;
%     velFlow(velFlow<velThresh)=0;
    
%     bwVelMap = velFlow;
%     bwVelMap(bwVelMap>velThresh) = 1;

%     bwVelMap=imbinarize(bwVelMap);
%     bwVelMap_filtered = bwareaopen(bwVelMap,60000);
%     streamPlotData = bwVelMap_filtered.*velocity_mag_val;
%     streamPlotData(:,:,10:end)=NaN;
    


%     [faces,verts] = isosurface(x_Rgrid3D,y_Rgrid3D,z_Rgrid3D,streamPlotData,0);
%     patch(ax2,'Faces',faces,'Vertices',verts,'FaceColor',[255, 0, 0]./255,'EdgeColor','None');
%     
    
    set(ax2,'ZDir','reverse');
    hold(ax2,'on');

    
    %% Eddy structure
    for differentEddyIndex=1:1:size(pathIndex,1)
        

        % Read eddy history
        eddyHistory = eddyPathHistory{pathIndex(differentEddyIndex)};
        eddyHistory=sortrows(eddyHistory,15);

        subpath=eddyPath{pathIndex(differentEddyIndex)};
        subGraphMap=subgraph(G,subpath);    
        firstEddyNode = int2str(eddyHistory(1,15))+"."+int2str(eddyHistory(1,16)+1);
        firstEddyNodeIndex=find(table2cell(subGraphMap.Nodes)==firstEddyNode);
        

        
        if(eddyHistory(1,15)>i)
            continue;
        else
            %% Current paths
            if(eddyHistory(end,15)>=i)   
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
            
            
            
            
            minX=min(eddyHistory(:,1));
            maxX=max(eddyHistory(:,1));
            minY=min(eddyHistory(:,2));
            maxY=max(eddyHistory(:,2));
            
            
            for index=1:1:eddyInThisFrameNum
                % Get each eddy in this eddy path history at this frame 
                % Plot the boundary of this eddy in 2D
                thisEddy = eddyInThisFrame(index,:);
                r = thisEddy(13)*resolution;  
                pos = [thisEddy(1),thisEddy(2)]; 
                t=0:0.001:(2*pi);  
                t=[t,0];
                eddy2DBoundary = plot(ax1,pos(1)+r*sin(t),pos(2)+r*cos(t),'k','linewidth',6);
                hold(ax1,'on');



                
                eddyFrameIndex = (eddyHistory(:,15)==i);
                
                if(eddyHistory(eddyFrameIndex,14)==1)
                    data = load(dataFilePath+"Seperated Structures/clockwise/Frame_"+num2str(FrameNum)+"_eddy_"+num2str(eddyHistory(eddyFrameIndex,16))+"_statistic.uocd");
                elseif(eddyHistory(eddyFrameIndex,14)==0)
                    data = load(dataFilePath+"Seperated Structures/counterclockwise/Frame_"+num2str(FrameNum)+"_eddy_"+num2str(eddyHistory(eddyFrameIndex,16))+"_statistic.uocd");
                else
                    error('Error: Can not find corresponding eddy data');
                end
                x = data(:,3);
                y = data(:,4);
                z = data(:,5);
                ow = data(:,6);
                u_obj = data(:,7);
                v_obj = data(:,8);
                temp = data(:,10);
                salt = data(:,11);
                
                % Get the boundary of this eddy structure
                bound = boundary(x,y,z,0.9);
                eddy3DBoundary = trisurf(bound,x,y,z,temp,'EdgeColor','none', 'FaceAlpha', '0.9','Parent', ax2);                
%                 localOb = trisurf(bound,x,y,z,sqrt(u_obj.^2+v_obj.^2),'EdgeColor','none', 'FaceAlpha', '0.9','Parent', ax2);
                localOb.SpecularExponent = 200;
                localOb.AmbientStrength = 0.8;
                hold(ax2,'on');
                

        
                if(i>1)
                    lastEddyNode = int2str(eddyHistory(eddyFrameIndex,15))+"."+int2str(eddyHistory(eddyFrameIndex,16)+1);
                    lastEddyNodeIndex=find(table2cell(subGraphMap.Nodes)==lastEddyNode);
                    thisEddyHistory=allpaths(subGraphMap,firstEddyNodeIndex,lastEddyNodeIndex);
                    thisEddyPath=eddyHistory(thisEddyHistory,[1 2 5 15]);
                    thisEddyPath = sortrows(thisEddyPath,4);


                    % Draw the trace of the eddy in 2D video
                    for historyIndex=1:1:size(thisEddyPath,1)-1
                        alphaValue=1-(size(thisEddyPath,1)-1-thisEddyPath(historyIndex,4)+thisEddyPath(1,4))*0.1;
                        if(alphaValue<0.4)
                            alphaValue=0.4;
                        end
                        eddy2DPath = plot(ax1,[thisEddyPath(historyIndex,1),thisEddyPath(historyIndex+1,1)], [thisEddyPath(historyIndex,2),thisEddyPath(historyIndex+1,2)],'Color',[0.4940 0.1840 0.5560,alphaValue],'linewidth',2);
            %                 plot(ax1,[CenterHistory(historyIndex,1),CenterHistory(historyIndex+1,1)], [CenterHistory(historyIndex,2),CenterHistory(historyIndex+1,2)],'Color',[0,0,0],'linewidth',2);

                        hold(ax1, 'on');
                    end
                    scatter(ax1,thisEddyPath(:,1), thisEddyPath(:,2),18,[0.4940 0.1840 0.5560],'filled');


                    %---------------------------------------------------------
                    % Add the code to draw both current and previous boundary
                    % (and arrows)
                    %---------------------------------------------------------                        

                    hold(ax1, 'on');

                    % Draw the trace of the eddy in 3D video
                    for historyIndex=1:1:size(thisEddyPath,1)-1
                        alphaValue=1-(size(thisEddyPath,1)-1-thisEddyPath(historyIndex,4)+thisEddyPath(1,4))*0.1;
                        if(alphaValue<0.4)
                            alphaValue=0.4;
                        end
                        eddy3DPath = plot3(ax2,[thisEddyPath(historyIndex,1),thisEddyPath(historyIndex+1,1)], [thisEddyPath(historyIndex,2),thisEddyPath(historyIndex+1,2)],[thisEddyPath(historyIndex,3),thisEddyPath(historyIndex+1,3)],'Color',[0.4940 0.1840 0.5560,alphaValue],'linewidth',4);
            %                 plot3(ax2,[CenterHistory(historyIndex,1),CenterHistory(historyIndex+1,1)], [CenterHistory(historyIndex,2),CenterHistory(historyIndex+1,2)],[CenterHistory(historyIndex,3),CenterHistory(historyIndex+1,3)],'Color',[0,0,0],'linewidth',2);
                        hold(ax2, 'on');
                    end
                    scatter3(ax2,thisEddyPath(:,1), thisEddyPath(:,2),thisEddyPath(:,3),8,[0.4940 0.1840 0.5560],'filled');
                    hold(ax2,'on');               
                end    
            end
    

            %% Past paths
%             else
%                 if(i>1)
%                     thisEddyNode=int2str(eddyHistory(end,15))+"."+int2str(eddyHistory(end,16)+1);
%                     thisEddyNodeIndex=find(table2cell(subGraphMap.Nodes)==thisEddyNode);
%                     thisEddyHistory=allpaths(subGraphMap,firstEddyNodeIndex,thisEddyNodeIndex);
%                     thisEddyPath=eddyHistory(thisEddyHistory,[1 2 5 15]);
%                     thisEddyPath = sortrows(thisEddyPath,4);
% 
% 
%                     % Draw the trace of the eddy in 2D video
%                     for historyIndex=1:1:size(thisEddyPath,1)-1
%                         alphaValue=1-(size(thisEddyPath,1)-1-thisEddyPath(historyIndex,4)+thisEddyPath(1,4))*0.1;
%                         if(alphaValue<0.4)
%                             alphaValue=0.4;
%                         end
%                         plot(ax1,[thisEddyPath(historyIndex,1),thisEddyPath(historyIndex+1,1)], [thisEddyPath(historyIndex,2),thisEddyPath(historyIndex+1,2)],'Color',[0.4940 0.1840 0.5560,0.3],'linewidth',2);
%             %                 plot(ax1,[CenterHistory(historyIndex,1),CenterHistory(historyIndex+1,1)], [CenterHistory(historyIndex,2),CenterHistory(historyIndex+1,2)],'Color',[0,0,0],'linewidth',2);
% 
%                         hold(ax1, 'on');
%                     end
%                     scatter(ax1,thisEddyPath(:,1), thisEddyPath(:,2),18,[0.4940 0.1840 0.5560],'filled');
%                     %---------------------------------------------------------
%                     % Add the code to draw both current and previous boundary
%                     % (and arrows)
%                     %---------------------------------------------------------       
%                     hold(ax1, 'on');
% 
%                     % Draw the trace of the eddy in 3D video
%                     for historyIndex=1:1:size(thisEddyPath,1)-1
%                         alphaValue=1-(size(thisEddyPath,1)-1-thisEddyPath(historyIndex,4)+thisEddyPath(1,4))*0.1;
%                         if(alphaValue<0.4)
%                             alphaValue=0.4;
%                         end
%                         plot3(ax2,[thisEddyPath(historyIndex,1),thisEddyPath(historyIndex+1,1)], [thisEddyPath(historyIndex,2),thisEddyPath(historyIndex+1,2)],[thisEddyPath(historyIndex,3),thisEddyPath(historyIndex+1,3)],'Color',[0.4940 0.1840 0.5560,0.3],'linewidth',4);
%             %                 plot3(ax2,[CenterHistory(historyIndex,1),CenterHistory(historyIndex+1,1)], [CenterHistory(historyIndex,2),CenterHistory(historyIndex+1,2)],[CenterHistory(historyIndex,3),CenterHistory(historyIndex+1,3)],'Color',[0,0,0],'linewidth',2);
%                         hold(ax2, 'on');
%                     end
%                     scatter3(ax2,thisEddyPath(:,1), thisEddyPath(:,2),thisEddyPath(:,3),8,[0.4940 0.1840 0.55600],'filled');
%                     hold(ax2,'on');               
%                 end
            end      
        end
%         for boundaryIndex = 1:1:size(B,1)
%             coastboundary=B{boundaryIndex};
%             plot(ax2,x_val(coastboundary(:,1)), y_val(coastboundary(:,2)),'k','LineWidth',4);
%             hold(ax2,'on');
%         end
    end
    

%     lightangle(ax2,-37.5,-45);


%     xlim(ax1,[minX-maxRadius-0.5,maxX+maxRadius+0.5]);
%     ylim(ax1,[minY-maxRadius-0.5,maxY+maxRadius+0.5]);

%     xlim(ax1,[-100,-50]);
%     ylim(ax1,[0,50]);

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


    coastTemp=temp_val(:,:,1);
    coastTemp(isnan(coastTemp))=100;
    coastTemp(~logical(coastTemp))=100;
    [coastPoints,coast]=contourf(ax2,x_Rgrid2D,y_Rgrid2D,coastTemp,[100 100],'DisplayName', 'land');
    hold(ax2,'on');    
    coast.FaceColor=[128, 85, 0]./255;

    cb1 = colorbar(ax1);
    cb1.Label.String="Temperature";
    cb.FontSize = 48;


    cb=colorbar(ax2);
%     colormap(ax2,"hot");
    cb.Label.String="Temperature";
    caxis(ax2,[34.5,38]);
    cb.FontSize = 48;
    colormap(ax2,jet);
 
    xlim(ax1,[-98,-85]);
    ylim(ax1,[18,30]);
    xlim(ax2,[-98,-85]);
    ylim(ax2,[18,30]);


%     zlim(ax2,[0,3300]);
%     camlight(ax2);
    lighting(ax2, 'gouraud');
    set(ax2,'ZDir','reverse');
    daspect(ax2,[1,1,800]);
    if(size(pathIndex)==1)
        title(ax2,"eddy "+int2str(pathIndex)+" from NA dataset Timestep "+num2str(FrameNum));
    else
        title(ax2,"eddies from NA dataset of Timestep "+num2str(FrameNum));
    end
    x2 = xlabel(ax2,'Longitude');
    y2 = ylabel(ax2,'Latitude');
    z2 = zlabel(ax2,'Depth');
    view(ax2,[48.34,25.78]);
    x2.FontSize = 48;
    y2.FontSize = 48;
    z2.FontSize = 48;
   

%     annotation(fh2,'textbox', [0.01 0.3 0.15 0.1], ...
%     'String', {'Yellow to blue structures:The detected eddies found by our proposed hybrid eddy detection algorithm.';' ';
%             'The blue-grey structure: seafloor';' ';
%             'Brown plane: land';' ';'Red isosurfaces: major ocean currents'}, ...
%     'Color', 'k', ...

%     'FontWeight', 'bold', ...
%     'FontSize',18, ...
%     'EdgeColor', 'none');
    
    if(i>1)
        legend(ax1,[eddy2DBoundary, eddy2DPath], "eddy 2D boundary", "eddy 2D path");
        legend(ax2,[eddy3DPath], "eddy 3D path");
    else
        legend(ax1,[eddy2DBoundary], "eddy 2D boundary");
%         legend(ax2,[eddy3DBoundary], "eddy 3D boundary");
    end

    F1=getframe(fh1);
    writeVideo(v,F1);
    F2=getframe(fh2);
    writeVideo(v2,F2);
    

end
close(v);
close(v2);
end

