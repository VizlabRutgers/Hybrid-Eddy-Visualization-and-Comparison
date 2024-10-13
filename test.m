
% clc;
% clear all;
% close all;

% creat the figure

fh1 = figure();
fh1.WindowState = 'maximized';
ax1=axes(fh1);

fh2 = figure();
fh2.WindowState = 'maximized';
ax2=axes(fh2);

% plot the seabed
% seabed = patch(ax2,'Faces',coastFaces, 'Vertices',coastVerts, ...  
% 'FaceColor',[.7,.7,.7],...
% 'EdgeColor','none',...
% 'FaceAlpha',1);
% seabed.AmbientStrength = 0.1;
% hold(ax2,'on');


srcData="../Dataset/NA_dataset/19980701.ocean_5day.nc";
x_val = double(ncread(srcData, "xh"));
y_val = double(ncread(srcData, "yh"));
z_val = double(ncread(srcData, "z_l"));

FrameNum=10;
startLoc = [1,1,1,FrameNum];
count = [length(x_val),length(y_val),length(z_val),1];
stride = [1,1,1,1];

startLoc_2D = [1,1,FrameNum];
count_2D = [length(x_val),length(y_val),1];
stride_2D = [1,1,1];

u_val = ncread(srcData, "u", startLoc, count, stride);
v_val = ncread(srcData, "v",startLoc, count, stride);
velMag_val=sqrt(u_val.^2+v_val.^2);
%     w_val = ncread(srcData, "W",startLoc, count, stride);
%     eta_val = ncread(srcData,"ETA",startLoc_2D, count_2D, stride_2D);
temp_val = ncread(srcData,"temp",startLoc, count, stride);
%     salinity_val = ncread(srcData,"salt",startLoc, count, stride);
[x_Rgrid2D,y_Rgrid2D] = ndgrid(x_val,y_val);
[x_Rgrid3D,y_Rgrid3D,z_Rgrid3D] = ndgrid(x_val,y_val,z_val);
meanVel = mean(mean(velMag_val(:,:,1),'omitnan'),'omitnan');

mean(mean(velMag_val));

% testVel = velMag_val(413,:,:);
testVel = velMag_val(413,:,:);
testVel = squeeze(testVel);
imagesc(testVel');
colorbar();

velThresh=0.2;
velFlow = velMag_val-mean(mean(velMag_val,'omitnan'),'omitnan');
velFlow(velFlow<velThresh)=0;

bwVelMap = velFlow;
bwVelMap(bwVelMap>velThresh) = 1;
bwVelMap=imbinarize(bwVelMap);
bwVelMap_filtered = bwareaopen(bwVelMap,3000);

f1=figure;
ax1=axes(f1);
[faces,verts] = isosurface(x_Rgrid3D,y_Rgrid3D,z_Rgrid3D,bwVelMap_filtered.*velMag_val,0);
patch(ax1,'Faces',faces,'Vertices',verts,'FaceColor','r','EdgeColor','None');
set(ax1,'ZDir','reverse');
view(3);

lighting gouraud;

% quiver(ax1,x_Rgrid2D,y_Rgrid2D,u_val(:,:,1),v_val(:,:,1),'r','AutoScaleFactor',3);
coastTemp=temp_val(:,:,1);
coastTemp(isnan(coastTemp))=100;
[coastPoints,coast]=contourf(ax2,x_Rgrid2D,y_Rgrid2D,coastTemp,[100 100]);

% other settings
set(ax2,'ZDir','reverse');

set(ax2,'YDir','normal');
cb=colorbar(ax2);
cb.Label.String="Salinity";
%         caxis(ax2,[-0.002,0]);
xlim(ax2,[-85,-65]);
ylim(ax2,[20,40]);
zlim(ax2,[0,1500]);
daspect(ax2,[1,1,300]);
title(ax2,"eddy from NA dataset");
ylabel(ax2,'Longitude');
xlabel(ax2,'Latitude');
zlabel(ax2,'Depth');
view(ax2,[35,35]);

% Read seabed data and object data
coastFaces = load('coastFaceData.mat','coastFaces').coastFaces;
coastVerts = load('coastVertData.mat','coastVerts').coastVerts;
bound=load('objBoundary.mat','bound').bound;

% set the light
camlight(ax2);
lighting(ax2, 'flat');


objIndex=["5.20","6.24","6.26"];
% objIndex=["5.20","6.26"];

for i =1:1:length(objIndex)
    color=[0,0,0];
    color(i)=1;
    objInfo=strsplit(objIndex(i),'.');
    objFrame=str2double(objInfo(1));
    objIndexInFrame=str2double(objInfo(2));
    historyCounter=find(allEddy(:,15)==objFrame & allEddy(:,16)==objIndexInFrame);
    
    if(allEddy(historyCounter,14)==1)
        data = load(dataFilePath+"Seperated Structures/clockwise/Frame_"+num2str(objFrame)+"_eddy_"+num2str(objIndexInFrame-1)+"_statistic.uocd");
    elseif(allEddy(historyCounter,14)==0)
        data = load(dataFilePath+"Seperated Structures/counterclockwise/Frame_"+num2str(objFrame)+"_eddy_"+num2str(objIndexInFrame-1)+"_statistic.uocd");
    else
        error('Error: Can not find corresponding eddy data');
    end
    x = data(:,3);
    y = data(:,4);
    z = data(:,5);
    salt = data(:,11);

    % plot the object
    bound = boundary(x,y,z,0.9);
    localOb = trisurf(bound,x,y,z,salt,'EdgeColor','none', 'FaceColor',color,'FaceAlpha', '0.9','Parent', ax2);
    localOb.SpecularExponent = 200;
    localOb.AmbientStrength = 0.8;
    hold(ax2,'on');

end








% uimage(ax2,x_val,y_val,temp_val(:,:,1)',"CDataMapping","scaled");

hold(ax2,'on');    

coast.FaceColor=[0.9290 0.6940 0.1250];








% Content below not validated
% Writing obj into trisurf
%----------------------------------------
% objTri = triangulation(bound,x,y,z);
% stlwrite(objTri, 'obj.stl','text');


%----------------------------------------