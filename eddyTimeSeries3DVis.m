function [] = eddyTimeSeries3DVis(srcData, property, eddyIndex, eddyPathHistory,dataFilePath)
%EDDYINDIVIDUAL3DVIS_INNER 
% data is the uocd style
% an axes is necessary for the visualization
% z_val is used to extract surface information

close all;

fh1 = figure();
ax1=axes(fh1);

view(3)
set(gca, "ZDir", "reverse");

colorbar();

x_val = ncread(srcData, property.x);
y_val = ncread(srcData, property.y);
z_val = ncread(srcData, property.z);



load("coastFaceData_NA.mat");
load("coastVertData_NA.mat");

% plot the seabed
% seabed = patch(ax1,'Faces',coastFaces, 'Vertices',coastVerts, ...  
% 'FaceColor',[.0,.0,.9],...
% 'EdgeColor','none',...
% 'FaceAlpha',0.1);

seabed = patch(ax1,'Faces',coastFaces, 'Vertices',coastVerts, ...  
'FaceColor',[.5,.5,.5],...
'EdgeColor','none',...
'FaceAlpha',0.3);
seabed.AmbientStrength = 0.1;
hold(ax1,'on');

lighting flat;



for differentEddyIndex=1:1:length(eddyIndex)
    

    % Read eddy history
    eddyHistory = eddyPathHistory{eddyIndex(differentEddyIndex)};
    eddyHistory=sortrows(eddyHistory,15);

    for index=[2,4,6,10,12,14]
    
            if(eddyHistory(index,14)==1)
                data = load(dataFilePath+"Seperated Structures/clockwise/Frame_"+num2str(eddyHistory(index,15))+"_eddy_"+num2str(eddyHistory(index,16))+"_statistic.uocd");
            elseif(eddyHistory(index,14)==0)
                data = load(dataFilePath+"Seperated Structures/counterclockwise/Frame_"+num2str(eddyHistory(index,15))+"_eddy_"+num2str(eddyHistory(index,16))+"_statistic.uocd");
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
    
%     
            surfaceData = data(round(data(:,5),4) == round(z_val(1),4),:);
            x_surface = surfaceData(:,3);
            y_surface = surfaceData(:,4);
%             
%             
%             
%             
%             
%             
            centroid_x = mean(x_surface);
            centroid_y = mean(y_surface);
            stretchMode = 1;
            
            if(stretchMode == 1)
                interpolatedSurfaceOffset_X = centroid_x*0.2;
                interpolatedSurfaceOffset_Y = centroid_y*0.2;    
            else
                interpolatedSurfaceOffset_X = index*3;
                interpolatedSurfaceOffset_Y = index*3;   
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
            opacity = index/size(eddyHistory,1);
            if(opacity<0.2) 
                opacity=0.2;
            end
            localOb = trisurf(bound,x,y,z,sqrt(u_obj.^2+v_obj.^2), 'FaceAlpha',opacity,'EdgeColor', 'none', 'Parent',ax1);
%             localOb = trimesh(bound,x,y,z,'facecolor'  'Parent',ax1);
            surface_boundary = boundary(x_surface,y_surface);
           hold on 
%            plot(x_surface(surface_boundary),y_surface(surface_boundary),'Color','#99D200','LineWidth',3);
           plot(x_surface(surface_boundary),y_surface(surface_boundary),'Color','k','LineWidth',3);

           hold on
    end
end

u_val = ncread(srcData,property.u,[1,1,1,eddyHistory(index,15)],[length(x_val), length(y_val),1,1],[1,1,1,1]);
v_val = ncread(srcData,property.v,[1,1,1,eddyHistory(index,15)],[length(x_val), length(y_val),1,1],[1,1,1,1]);

hold on;
% quiver(x_val,y_val,u_val', v_val');

[x_Rgrid2D,y_Rgrid2D] = meshgrid(x_val,y_val);
[x_start,y_start] = meshgrid(x_val(1):1:x_val(end),y_val(1):1:y_val(end));
hold on
% streamline(x_Rgrid2D, y_Rgrid2D,u_val', v_val', x_start,y_start);
hold on

xlim([-95,-80]);
ylim([23,30]);
zlim([0,3000]);
daspect([1,1,500]);
xlabel("longitude");
ylabel("latitude");
zlabel("depth");
colormap("hot");
cb = colorbar();
cb.Label.String = "velocity magnitude";
set(gca,"FontSize",24);