function [] = eddyIndividual3DVis_inner(data,ax2,z_val,stretchMode,index,srcData, property,horizontalOffset, background_option)
%EDDYINDIVIDUAL3DVIS_INNER 
% data is the uocd style
% an axes is necessary for the visualization
% z_val is used to extract surface information



x = data(:,3);
y = data(:,4);
z = data(:,5);
ow = data(:,6);
u_obj = data(:,7);
v_obj = data(:,8);
temp = data(:,10);
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
y_lessRange = minY;



centroid_x = mean(x_surface);
centroid_y = mean(y_surface);


if(stretchMode == 1)
    interpolatedSurfaceOffset_X = centroid_x*9;
    interpolatedSurfaceOffset_Y = centroid_y*9;    
else
    interpolatedSurfaceOffset_X = (index-1)*6-centroid_x;
    interpolatedSurfaceOffset_Y = -centroid_y + horizontalOffset;   
end

x = x+interpolatedSurfaceOffset_X;
y = y+interpolatedSurfaceOffset_Y;



% Get the boundary of this eddy structure
bound = boundary(x,y,z,0.9);
% edgeAlpha = index / wholePathSize;
% edgeAlpha = edgeAlpha*0.7;
% if(edgeAlpha<0.3) 
%     edgeAlpha=0.3;
% end

x_val = ncread(srcData, property.x);
y_val = ncread(srcData,property.y);
temp = ncread(srcData, property.temp,[1,1,1,1], [500,500,1,1]);
temp =temp';
[xGrid, yGrid] = meshgrid(x_val,y_val);
mask = (xGrid < x_lessRange | xGrid > x_overRange) | (yGrid < y_lessRange | yGrid > y_overRange);

temp(mask) = nan;



img=zeros([size(temp),3]);
img_r = img(:,:,1);
img_g = img(:,:,2);
img_b = img(:,:,3);
img_r(temp==0) = 0.55;
img_g(temp==0) = 0.45;
img_b(temp==0) = 0.25;
img(:,:,1) = img_r;
img(:,:,2) = img_g;
img(:,:,3) = img_b;
alphaData = temp==0;



edgeAlpha=0;
localOb = trisurf(bound,x,y,z,sqrt(u_obj.^2+v_obj.^2), 'FaceAlpha', '0.8', 'EdgeAlpha', edgeAlpha,'Parent',ax2);


hold on

if(background_option == 1)
    
    land=imshow(img,'XData',x_val+ interpolatedSurfaceOffset_X,'YData',y_val+ interpolatedSurfaceOffset_Y);
    set(land, 'AlphaData', alphaData);
    axis on
    hold on
    
    coastFaces = load('coastFaceData_RS.mat','coastFaces').coastFaces;
    coastVerts = load('coastVertData_RS.mat','coastVerts').coastVerts;
    vertX = coastVerts(:,1);
    vertY = coastVerts(:,2);

    x_faces = vertX(coastFaces);
    y_faces = vertY(coastFaces);

    in_x_range = all(x_faces > x_lessRange & x_faces < x_overRange, 2);
    in_y_range = all(y_faces > y_lessRange & y_faces < y_overRange, 2);

    face_mask = in_x_range & in_y_range;
    coastFaces_sub = coastFaces(face_mask, :);


    coastVerts_offset = coastVerts;
    coastVerts_offset(:,1) = coastVerts_offset(:,1) + interpolatedSurfaceOffset_X;
    coastVerts_offset(:,2) = coastVerts_offset(:,2) + interpolatedSurfaceOffset_Y;



    seabed = patch(ax2,'Faces',coastFaces_sub, 'Vertices',coastVerts_offset, ...  
    'FaceColor',[0.7, 0.7, 0.7],...
    'EdgeColor','none',...
    'FaceAlpha',0.85);

    hold on
    h_bg = patch(NaN, NaN, [0.55,0.45,0.25],'EdgeColor','none');

    legend([seabed,h_bg],{"Seabed","sand"});
end
view(3);


end

