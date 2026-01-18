function [] = eddyBioVis_inner(data,ax2,z_val,stretchMode,index,srcData, property)
%EDDYINDIVIDUAL3DVIS_INNER 
% data is the uocd style
% an axes is necessary for the visualization
% z_val is used to extract surface information

coastFaces = load('coastFaceData_RS.mat','coastFaces').coastFaces;
coastVerts = load('coastVertData_RS.mat','coastVerts').coastVerts;

x = data(:,3);
y = data(:,4);
z = data(:,5);
ow = data(:,6);
u_obj = data(:,7);
v_obj = data(:,8);
temp_withRegionMask = data(:,10);
salt = data(:,11);


surfaceData = data(round(data(:,5),4) == round(z_val(1),4),:);
x_surface = surfaceData(:,3);
y_surface = surfaceData(:,4);

minX = min(x);
maxX = max(x);
minY = min(y);
maxY = max(y);

x_overRange = maxX+0.5;
x_lessRange = minX-2;
y_overRange = maxY+2.5;
y_lessRange = minY+0.5;


regionX_min = min(x_lessRange, minX-1);
regionX_max = max(x_lessRange, maxX+1);
regionY_min = min(y_lessRange, minY-1);
regionY_max = max(y_lessRange, maxY+1);


centroid_x = mean(x_surface);
centroid_y = mean(y_surface);


if(stretchMode == 1)
    interpolatedSurfaceOffset_X = centroid_x*9;
    interpolatedSurfaceOffset_Y = centroid_y*9;    
else
    interpolatedSurfaceOffset_X = (index-1)*6-centroid_x;
    interpolatedSurfaceOffset_Y = -centroid_y;   
end

% x = x+interpolatedSurfaceOffset_X;
% y = y+interpolatedSurfaceOffset_Y;



% Get the boundary of this eddy structure
% bound = boundary(x,y,z,0.9);
% edgeAlpha = index / wholePathSize;
% edgeAlpha = edgeAlpha*0.7;
% if(edgeAlpha<0.3) 
%     edgeAlpha=0.3;
% end

x_val = ncread(srcData, property.x);
y_val = ncread(srcData,property.y);
temp = ncread(srcData, property.temp,[1,1,1,1], [length(x_val), length(y_val),1,1]);
temp =temp';
temp_withRegionMask = temp;
[xGrid, yGrid] = meshgrid(x_val,y_val);
mask = (xGrid < x_lessRange | xGrid > x_overRange) | (yGrid < y_lessRange | yGrid > y_overRange);

oceansurfaceMask = (xGrid < regionX_min | xGrid > regionX_max) | (yGrid < regionY_min | yGrid > regionY_max);
temp_withRegionMask(mask) = nan;



img=zeros([size(temp_withRegionMask),3]);
img_r = img(:,:,1);
img_g = img(:,:,2);
img_b = img(:,:,3);
img_r(temp_withRegionMask==0) = 0.55;
img_g(temp_withRegionMask==0) = 0.45;
img_b(temp_withRegionMask==0) = 0.25;
img(:,:,1) = img_r;
img(:,:,2) = img_g;
img(:,:,3) = img_b;
alphaData = temp_withRegionMask==0;

land=imshow(img,'XData',x_val,'YData',y_val,'Parent',ax2);
set(land, 'AlphaData', alphaData);
axis on
hold on

edgeAlpha=0;
% localOb = trisurf(bound,x,y,z,sqrt(u_obj.^2+v_obj.^2), 'FaceAlpha', '0.8', 'EdgeAlpha', edgeAlpha,'Parent',ax2);

temp_withOceanMask = temp;
temp_withOceanMask(oceansurfaceMask) = nan;
temp_withOceanMask(temp_withOceanMask==0) = nan;
% temp_surface = imagesc(ax2,x_val, y_val, temp_withOceanMask);


alpha = 0.5 * ones(size(temp_withOceanMask));   % 正常区域设为透明度 0.1
alpha(isnan(temp_withOceanMask)) = 0;           % NaN 区域完全透明

% set(temp_surface, 'AlphaData', alpha);  % 让 NaN 区域透明
hold(ax2, "on")


centerX = mean(surfaceData(:,3));
centerY = mean(surfaceData(:,4));

pos = [centerX,centerY]; 
t=0:0.001:(2*pi);  
t=[t,0];
radius =(max(surfaceData(:,3)) - min(surfaceData(:,3)))/2;  
depth = ones(1,length(t))*z_val(1);
interpolation_BoundaryPlot_handle = plot3(ax2, pos(1)+radius*sin(t),pos(2)+radius*cos(t), depth,'r','linewidth',3);
hold(ax2, "on")

vertX = coastVerts(:,1);
vertY = coastVerts(:,2);

x_faces = vertX(coastFaces);
y_faces = vertY(coastFaces);

in_x_range = all(x_faces > x_lessRange & x_faces < x_overRange, 2);
in_y_range = all(y_faces > y_lessRange & y_faces < y_overRange, 2);

face_mask = in_x_range & in_y_range;
coastFaces_sub = coastFaces(face_mask, :);


coastVerts_offset = coastVerts;
coastVerts_offset(:,1) = coastVerts_offset(:,1);
coastVerts_offset(:,2) = coastVerts_offset(:,2);



seabed = patch(ax2,'Faces',coastFaces_sub, 'Vertices',coastVerts_offset, ...  
'FaceColor',[0.7, 0.7, 0.7],...
'EdgeColor','none',...
'FaceAlpha',0.85);

hold(ax2, "on");
h_bg = patch(NaN, NaN, [0.55,0.45,0.25],'EdgeColor','none');

% legend(ax2, [seabed,h_bg,interpolation_BoundaryPlot_handle],{"Seabed","land","eddy boundary on the surface"});
view(ax2,3);

xlim(ax2,[regionX_min,regionX_max]);
ylim(ax2,[regionY_min,regionY_max]);
end

