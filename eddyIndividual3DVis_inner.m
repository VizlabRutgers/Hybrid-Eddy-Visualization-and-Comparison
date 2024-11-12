function [] = eddyIndividual3DVis_inner(data,ax2,z_val,stretchMode,index,wholePathSize)
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






        centroid_x = mean(x_surface);
        centroid_y = mean(y_surface);


        if(stretchMode == 1)
            interpolatedSurfaceOffset_X = centroid_x*9;
            interpolatedSurfaceOffset_Y = centroid_y*9;    
        else
            interpolatedSurfaceOffset_X = index*3;
            interpolatedSurfaceOffset_Y = index*3;   
        end

        x = x+interpolatedSurfaceOffset_X;
        y = y+interpolatedSurfaceOffset_Y;

        
        
        % Get the boundary of this eddy structure
        bound = boundary(x,y,z,0.9);
        edgeAlpha = index / wholePathSize;
        edgeAlpha = edgeAlpha*0.7;
        if(edgeAlpha<0.3) 
            edgeAlpha=0.3;
        end
        localOb = trimesh(bound,x,y,z,sqrt(u_obj.^2+v_obj.^2), 'FaceAlpha', '0.5', 'EdgeAlpha', edgeAlpha,'Parent',ax2);

end

