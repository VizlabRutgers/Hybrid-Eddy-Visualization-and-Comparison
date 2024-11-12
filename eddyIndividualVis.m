function [outputArg1,outputArg2] = eddyIndividualVis(pathIndex,G,eddyPath,eddyPathHistory,dataFilePath,srcData, property, eddyIndex,stretchMode)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
close all;

% Video Recorder
% v = VideoWriter('eddyHistory_2D_redsea.mp4');
% v.FrameRate=1;
% v.Quality=100;
% open(v);
% 
% v2 = VideoWriter('eddyHistory_3D_redsea.mp4');
% v2.FrameRate=1;
% v2.Quality=100;
% open(v2);


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




% Compute seabed data

% coastRegion = ncread(srcData,property.temp,[1,1,1,1], [length(x_val),length(y_val),length(z_val),1], [1,1,1,1]);
% coast2D=~(coastRegion(:,:,1));
% [B,L] = bwboundaries(coast2D,'holes');
% coastRegion(~logical(coastRegion))=-1;
% [coastFaces,coastVerts]=isosurface(x_val,y_val,z_val,permute(coastRegion,[2,1,3]),-1);
% 
% save('coastFaceData_RS.mat','coastFaces');
% save('coastVertData_RS.mat','coastVerts');

% Load seabed data
% coastFaces = load('coastFaceData_NA.mat','coastFaces').coastFaces;
% coastVerts = load('coastVertData_NA.mat','coastVerts').coastVerts;

% coastFaces = load('coastFaceData_NP.mat','coastFaces').coastFaces;
% coastVerts = load('coastVertData_NP.mat','coastVerts').coastVerts;

% coastFaces = load('coastFaceData_RS.mat','coastFaces').coastFaces;
% coastVerts = load('coastVertData_RS.mat','coastVerts').coastVerts;

maxFrame = cellfun(@(x) size(x,1), eddyPathHistory);
maxFrame = max(maxFrame);

% Here the i indicates the time frame in eddy history

% for i = 1
%     if (i == maxFrame)
%         i = maxFrame;
%     end

FrameNum=1;
startLoc = [1,1,1,FrameNum];
count = [length(x_val),length(y_val),length(z_val),1];
stride = [1,1,1,1];

startLoc_2D = [1,1,FrameNum];
count_2D = [length(x_val),length(y_val),1];
stride_2D = [1,1,1];


% u_val = ncread(srcData, property.u, startLoc, count, stride);
% v_val = ncread(srcData, property.v,startLoc, count, stride);
%     w_val = ncread(srcData, "W",startLoc, count, stride);
%     eta_val = ncread(srcData,"ETA",startLoc_2D, count_2D, stride_2D);
% temp_val = ncread(srcData,property.temp,startLoc, count, stride);
%     salinity_val = ncread(srcData,"salt",startLoc, count, stride);
[x_Rgrid3D,y_Rgrid3D,z_Rgrid3D] = ndgrid(x_val,y_val,z_val);
% velocity_mag_val = (sqrt(u_val.^2 +v_val.^2));
% 
% 
% % Get quvier map
% [x_Rgrid2D,y_Rgrid2D] = ndgrid(x_val(startLoc(1):1:(startLoc(1)+count(1))-1),y_val(startLoc(2):1:(startLoc(2)+count(2))-1));
% % 
% quiver(ax1,x_Rgrid2D,y_Rgrid2D,u_val(:,:,1),v_val(:,:,1),'r','AutoScaleFactor',3);
% hold(ax1,'on');

% % Draw the sea bed
% seabed = patch(ax2,'Faces',coastFaces, 'Vertices',coastVerts, ...  
% 'FaceColor',[.8, .8, .8],...
% 'EdgeColor','none',...
% 'FaceAlpha',0.5);
% isonormals(temp_val, seabed);
% % seabed.SpecularColorReflectance = 0;
% %     seabed.SpecularExponent = 80;
% seabed.AmbientStrength = 0.1;





% 
% 
% for i=1:1:length(eddyHistory)
%     pathLength(i) = size(eddyHistory{i},1);
% end
% 
% [loc,val] = max(pathLength);

% velThresh=0.8;
% velFlow = (velocity_mag_val-mean(mean(velocity_mag_val,'omitnan'),'omitnan'))./velocity_mag_val;
% velFlow(velFlow<velThresh)=0;
% 
% bwVelMap = velFlow;
% bwVelMap(bwVelMap>velThresh) = 1;
% 
% bwVelMap=imbinarize(bwVelMap);
% bwVelMap_filtered = bwareaopen(bwVelMap,60000);
% streamPlotData = bwVelMap_filtered.*velocity_mag_val;
% streamPlotData(:,:,10:end)=NaN;
% 
% 
% [faces,verts] = isosurface(x_Rgrid3D,y_Rgrid3D,z_Rgrid3D,streamPlotData,0);
% patch(ax2,'Faces',faces,'Vertices',verts,'FaceColor',[255, 0, 0]./255,'EdgeColor','None');
% set(ax2,'ZDir','reverse');


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

% for differentEddyIndex=1:1:size(eddyPathHistory,1)
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
    
    
    wholePathSize = size(eddyHistory,1);
    for index=1:1:wholePathSize
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
        r = thisEddy(13)*resolution;  
        pos = [thisEddy(1),thisEddy(2)]; 
        t=0:0.001:(2*pi);  
        t=[t,0];
        BoundaryPlot_handle = plot(ax1,pos(1)+r*sin(t),pos(2)+r*cos(t),'k','linewidth',6);
        hold(ax1,'on');

        eddyIndividual3DVis_inner(data,ax2,z_val,stretchMode,index,wholePathSize);
    end      
end

set(ax2,'ZDir','reverse');
view(3)
daspect([1 1 100])


grid on;

% xticklabels("")
% yticklabels("")
% zticklabels("")


camlight(ax2);
lighting(ax2, 'flat');
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


cb=colorbar(ax2);
%     colormap(ax2,"hot");
cb.Label.String="Net Velocity";
caxis(ax2,[0,0.3]);
cb.FontSize = 48;

% xlim(ax2,[35,38]);
% ylim(ax2,[21,25]);
% zlim(ax2,[0,50]);
set(ax2,'ZDir','reverse');
% daspect(ax2,[1,1,800]);

% pbaspect(ax2,[1,1,1])

if(size(pathIndex)==1)
    title(ax2,"eddy "+int2str(pathIndex)+" from NA dataset Frame"+num2str(FrameNum));
else
    title(ax2,"eddies from NA dataset of Day "+num2str((FrameNum-1)*5+1));
end
x2 = xlabel(ax2,'Longitude');
y2 = ylabel(ax2,'Latitude');
z2 = zlabel(ax2,'Depth');
view(ax2,[68.2,24]);
x2.FontSize = 48;
y2.FontSize = 48;
z2.FontSize = 48;


annotation(fh2,'textbox', [0.01 0.3 0.15 0.1], ...
'String', {'Yellow to blue structures:The detected eddies found by our proposed hybrid eddy detection algorithm.';' ';
        'The blue-grey structure: seafloor';' ';
        'Brown plane: land';' ';'Red isosurfaces: major ocean currents'}, ...
'Color', 'k', ...
'FontWeight', 'bold', ...
'FontSize',18, ...
'EdgeColor', 'none');

%% costal sea bed
% temp_val = ncread(srcData,property.temp,startLoc, count, stride);
% coastTemp=temp_val(:,:,1);
% coastTemp(isnan(coastTemp))=100;
% coastTemp(~logical(coastTemp))=100;
% [coastPoints,coast]=contourf(ax2,x_Rgrid2D,y_Rgrid2D,coastTemp,[100 100]);
% hold(ax2,'on');    
% coast.FaceColor=[128, 85, 0]./255;




% F1=getframe(fh1);
% writeVideo(v,F1);
% F2=getframe(fh2);
% writeVideo(v2,F2);
%     
%     
% close(v);
% close(v2);


end

