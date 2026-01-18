function [outputArg1,outputArg2] = windingAngle_HorizontalComparison(x_val, y_val, x_min, y_min, x_max, y_max)
lon = x_val;
lat = y_val;
lon_max = x_max;
lon_min = x_min;
lat_max = y_max;
lat_min = y_min;

if plt_debug == 1 % ellipse plots
    [lats,lons] = meshgrid(lat,lon);
    lon_inLim = lon(lon<=lon_max & lon>=lon_min);
    lat_inLim = lat(lat<=lat_max & lon>=lat_min);

    if lonlatlimits == 1
        j = lons >= lon_min & lons <= lon_max & lats >= lat_min & lats <= lat_max;
        u = u(j);
        v = v(j);
        lons = lons(j);
        lats = lats(j);

        j = lon >= lon_min & lon <= lon_max;
        lon1 = lon(j);
        j = lat >= lat_min & lat <= lat_max;
        lat1 = lat(j);
        [lats1,lons1] = meshgrid(lat1,lon1);
        lons = reshape(lons,size(lats1));
        lats = reshape(lats,size(lats1));
        u = reshape(u,size(lats1));
        v = reshape(v,size(lats1));

    else
        lon_min = double(floor(min(lon)));
        lon_max = double(ceil(max(lon)));
        lat_min = double(floor(min(lat)));
        lat_max = double(ceil(max(lat)));
    end

        %             figure
    %             m_quiver(lons,lats,u,v)
    %             axis equal
    
    h = figure('color','w','renderer','painters');
    m_proj('Mercator','lon',[lon_min lon_max],'lat',[lat_min lat_max])
    % m_gshhs_h('patch',[0.8 0.8 0.8]); % coastline
    if lonlatlimits == 1 && highrescoast == 1
        m_gshhs_i('patch',[0.8 0.8 0.8]);
    else
        m_gshhs_l('patch',[0.8 0.8 0.8]);
    end
    hold on
    % levels = [-10:-5:-30 -40 -50 -100 -150 -200]; % isobaths
    % levels = [-100 -1000]; % isobaths
    % m_etopo2('contour',levels,'edgecolor','k','linewidth',1,'linestyle','--','ShowText','on');

%     hold on
    quiver_handle = m_quiver(lons,lats,u,v);
    hold on

%     m_pcolor(lons, lats, sqrt(u.^2+v.^2));


    
    hold on
    m_grid('box','fancy','tickdir','out','linewidth',1,'linestyle','none','fontsize',14);  % nice grid
%     cb = colorbar();
%     cb.Label.String = "Velocity Magnitude";
    
    hold on
    
    for eddyIndex = 1:1:eddies
        %     d(:,1);
        %     d(:,2);



        origin_center = m_plot(eddy_centerX_original(eddyIndex),eddy_centerY_original(eddyIndex),...
            'kx','linewidth',4,'markersize',12);
        hold on
        improved_center = m_plot(eddy_centerX_new(eddyIndex),eddy_centerY_new(eddyIndex),...
            'rx','linewidth',4,'markersize',12);
        hold on
%         m_plot(eddy_centerX_boundary(eddies) + ellipse_x_boundary,...
%             eddy_centerY_boundary(eddies) + ellipse_y_boundary,'r--','linewidth',2)
        origin_boundary = m_plot(eddy_centerX_original(eddyIndex) + ellipse_x_original{eddyIndex},...
            eddy_centerY_original(eddyIndex) + ellipse_y_original{eddyIndex},'k-','linewidth',2 );
        hold on
%         max_vel_streamline = m_plot(eddy_boundary_maxVel_x{eddyIndex},...
%             eddy_boundary_maxVel_y{eddyIndex},'r--','linewidth',2);
%         hold on
        improved_boundary = m_plot(eddy_boundary_maxVel_x_closed{eddyIndex},...
            eddy_boundary_maxVel_y_closed{eddyIndex},'r-','linewidth',2);
        hold on
    end

%     legend([origin_center origin_boundary],...,
%     {'original winding angle center', 'original winding angle boundary'}, 'Location', 'best');


    legend([origin_center improved_center origin_boundary improved_boundary, quiver_handle],...,
        {'original winding angle center', 'improved winding angle center', 'original winding angle boundary' ...
        , 'improved winding angle boundary', 'velocity field'}, 'Location', 'best');
    % show time of calculation
    daspect([1 1 1]);
    set(gca, "FontSize", 20);
    xlabel("Longitude");
    ylabel("Latitude");
    
end


if eddies == 0
    eddy_centerX_boundary = [];
    eddy_centerY_boundary = [];
    eddy_dir = [];
    eddy_angular_vel_boundary = [];
    eddy_streamlines = [];
    eddy_length_x_boundary = [];
    eddy_length_y_boundary = [];
    eddy_ellipse_theta_boundary = [];
end
% save data
save(snmat,'eddy_centerX_boundary','eddy_centerY_boundary',...
    'eddy_dir','eddy_angular_vel_boundary','eddy_streamlines',...
    'eddy_length_x_boundary','eddy_length_y_boundary','eddy_ellipse_theta_boundary');

if plt_debug == 1 % save plot
    savefig(sn)
    print(sn,'-dpng')
    if numfiles ~= 1
        close(h)
    end
end

disp([num2str(eddies) ' eddies found'])
for i_eddy = 1:eddies
    disp(['eddy #' num2str(i_eddy) ' has ' ...
        num2str(eddy_streamlines(i_eddy)) ' streamlines'])
end
disp(' ')

end

