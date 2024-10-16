function [] = eddyGlobalStat(allEddy, FrameLimit, srcData,sizeLimit,patchSize, property,FrameNum)
% allPath: All of eddy's path (used for containg every eddy here)
% FrameNum: Decide which frame to be counted in the global statistic
% srcData: The src data of .nc file
% sizeLimit: Filter out the eddy lower than this limit
% patchSize: Decide the size of patch
%   Detailed explanation goes here
%------------------------------------
    close all;

    % Read data
    x_val = double(ncread(srcData, property.x));
    y_val = double(ncread(srcData, property.y));
    z_val = double(ncread(srcData, property.z));
    startLoc = [1,1,1,1];
    count = [length(x_val),length(y_val),length(z_val),1];
    stride = [1,1,1,1];
    
    startLoc_2D = [1,1,1];
    count_2D = [length(x_val),length(y_val),1];
    stride_2D = [1,1,1];



    u_val = ncread(srcData, property.u, startLoc, count, stride);
    v_val = ncread(srcData, property.v,startLoc, count, stride);
%     w_val = ncread(srcData, property.w,startLoc, count, stride);
    eta_val = ncread(srcData,"ETA",startLoc_2D, count_2D, stride_2D);
    temp_val = ncread(srcData,property.temp,startLoc, count, stride);
%     salinity_val = ncread(srcData,"salt",startLoc, count, stride);

    xSize = length(x_val);
    ySize = length(y_val);
    
    if(isempty(FrameNum))
        FrameNum = FrameLimit.min:1:FrameLimit.max;
    end

    % Get the coast boundary
    coastRegion=~isnan(temp_val(:,:,1));
    [B,L] = bwboundaries(coastRegion,'holes');
    
    xmin = x_val(1);
    ymin = y_val(1);
    
    velocity_mag_val = (sqrt(u_val.^2 +v_val.^2));
    
    
    u_val_test = u_val(1:xSize-1,:,:) - u_val(2:xSize,:,:);
    u_x = u_val;
    u_x(1:xSize-1,:,:) = u_val_test;
    u_val_test = u_val(:,1:ySize-1,:) - u_val(:,2:ySize,:);
    u_y = u_val;
    u_y(:,1:ySize-1,:) = u_val_test;
    
    v_val_test = v_val(1:xSize-1,:,:) - v_val(2:xSize,:,:);
    v_x = v_val;
    v_x(1:xSize-1,:,:) = v_val_test;
    v_val_test = v_val(:,1:ySize-1,:) - v_val(:,2:ySize,:);
    v_y = v_val;
    v_y(:,1:ySize-1,:) = v_val_test;
    
    s_n = (u_x-v_y).^2;
    s_s = (u_y+v_x).^2;
    omega = (v_x-u_y).^2;
    
    s_n(s_n==0) = NaN;
    s_s(s_n==0) = NaN;
    omega(s_n==0) = NaN;
    ow_val = s_s+s_n-omega;
    ow_val = ow_val./(velocity_mag_val.^2);


    % Initialize the patch array foe the statistic data
    xSize = x_val(1:patchSize:end);
    ySize = y_val(1:patchSize:end);
%     [xGrid,yGrid]=ndgrid(xSize,ySize);
%     xRound = round(x_val,4);
%     yRound = round(y_val,4);
    eddyNumberDivision=zeros(length(xSize),length(ySize));
    eddySizeDivision=zeros(length(xSize),length(ySize));
    eddyOWDivision=zeros(length(xSize),length(ySize));
    eddyNumTimeSizeDivision=zeros(length(xSize),length(ySize));
    
    % Tranverse every frame
    for frameIndex=FrameNum
        eddyPath=allEddy(allEddy(:,15)==frameIndex & allEddy(:,13)>=sizeLimit.lower,:);
        
%         Code below is used to see statistic of each eddy without patch
%         setting
%         ----------------------------------------------------------
        figure,
        ax1 = axes;


%         uimage(ax1,x_val,y_val,~temp_val(:,:,1)',"CDataMapping","scaled");
%         color_map=[1 1 1; 0.3 0.3 0.3];
%         colormap(color_map);
        
        hold on
        [x_Rgrid2D,y_Rgrid2D] = ndgrid(x_val(startLoc(1):1:(startLoc(1)+count(1))-1),y_val(startLoc(2):1:(startLoc(2)+count(2))-1));
        % 
        quiver(ax1,x_Rgrid2D,y_Rgrid2D,u_val(:,:,1),v_val(:,:,1),10,'r','AutoScaleFactor',0.5);


        title("Eddy Size for frame"+num2str(1));
        ylabel('Latitude');
        xlabel('Longitude');
        daspect([1,1,1]);
        hold on
        scatter(ax1,eddyPath(:,1),eddyPath(:,2),[],eddyPath(:,13),"filled");
        cb=colorbar();
        cb.Label.String='Eddy Size';
        caxis([4,10]);
        set(ax1,'YDir','normal');

    
    
        figure,   
        ax2 = axes;
        uimage(ax2,x_val(startLoc(1):1:(startLoc(1)+count(1))-1),y_val(startLoc(2):1:(startLoc(2)+count(2))-1),eta_val',"CDataMapping","scaled");


        title("SSH for frame"+num2str(1));
        ylabel('Latitude');
        xlabel('Longitude');
        daspect([1,1,1]);
        hold on


        scatter(ax2,eddyPath(:,1),eddyPath(:,2),'k',"filled");
        cb=colorbar();
        cb.Label.String='SSH';
        set(ax2,'YDir','normal');


        figure,  
        ax3 = axes;
        uimage(ax3,x_val(startLoc(1):1:(startLoc(1)+count(1))-1),y_val(startLoc(2):1:(startLoc(2)+count(2))-1),omega(:,:,1)',"CDataMapping","scaled");


        title("Vorticity for frame"+num2str(1)+" on surface");
        ylabel('Latitude');
        xlabel('Longitude');
        daspect([1,1,1]);
        hold on


        scatter(ax3,eddyPath(:,1),eddyPath(:,2),'k',"filled");
        cb=colorbar();
        cb.Label.String='Relative Vorticity';
        set(ax3,'YDir','normal');
        caxis([0,0.01]);




        for i = 1:1:length(size(eddyPath,1))
            thisEddy=eddyPath(i,:);
            r = thisEddy(13)*0.04;  
            pos = [thisEddy(1),thisEddy(2)]; 
            t=0:0.001:(2*pi);  
            t=[t,0];
            BoundaryPlot_handle = plot(pos(1)+r*sin(t),pos(2)+r*cos(t));
            hold on
            legend([BoundaryPlot_handle(1)],"Eddy Boundary");


        end

%         patch(pos(1)+r*sin(t),pos(2)+r*cos(t), h*ones(size(t)),'k');
%         ----------------------------------------------------------




        % Clear the temp data of patch array in each loop(frame)
        eddySizeDivisionTemp=cell(length(xSize),length(ySize));
        eddySizeDivisionTemp=cellfun(@(x) [],eddySizeDivisionTemp,'UniformOutput',false);

        eddyNumberDivisionTemp=cell(length(xSize),length(ySize));
        eddyNumberDivisionTemp=cellfun(@(x) [],eddyNumberDivisionTemp,'UniformOutput',false);
        
        eddyOWDivisionTemp=cell(length(xSize),length(ySize));
        eddyOWDivisionTemp=cellfun(@(x) [],eddyOWDivisionTemp,'UniformOutput',false);        
        
        % Find out which patch this eddy belongs to. Take the location of
        % eddy minus the coordinate of the patch array and then find the first
        % positive result
        for i = 1:1:size(eddyPath,1)
            currentEddy=eddyPath(i,:);
            xLoc=xSize-currentEddy(1);
            xLoc=find(xLoc>0,1,'first');
            if(isempty(xLoc))
                xLoc=length(xSize);
            end
            yLoc=ySize-currentEddy(2);
            yLoc=find(yLoc>0,1,'first');
            if(isempty(yLoc))
                yLoc=length(ySize);
            end
            eddySizeDivisionTemp{xLoc,yLoc}=[eddySizeDivisionTemp{xLoc,yLoc},currentEddy(13)];
%             [d,eddyXLoc]=min(abs(currentEddy(1)-x_val));
%             [d,eddyYLoc]=min(abs(currentEddy(2)-y_val));
            eddyOWDivisionTemp{xLoc,yLoc}=[eddyOWDivisionTemp{xLoc,yLoc},currentEddy(9)];
        end
        
        % Arithmetic function to get mean/sum/size of different statistic
        % result
        eddyNumberDivisionTemp=cellfun(@(x) length(x),eddySizeDivisionTemp);
        eddyNumTimeSizeTemp=cellfun(@(x) sum(x),eddySizeDivisionTemp);
        eddySizeDivisionTemp=cellfun(@(x) mean(x),eddySizeDivisionTemp);
        eddyOWDivisionTemp=cellfun(@(x) mean(x),eddyOWDivisionTemp);  

        
        eddySizeDivisionTemp(isnan(eddySizeDivisionTemp))=0;
        eddyOWDivisionTemp(isnan(eddyOWDivisionTemp))=0;

        % Put the temp data in each loop(frame) to the overall data
        eddyNumberDivision=eddyNumberDivision+eddyNumberDivisionTemp;
        eddySizeDivision=eddySizeDivision+eddySizeDivisionTemp;
        eddyOWDivision=eddyOWDivision+eddyOWDivisionTemp;
        eddyNumTimeSizeDivision=eddyNumTimeSizeDivision+eddyNumTimeSizeTemp;
    
    end
    
%     eddyNumberDivision=cellfun(@(x) length(x),eddySizeDivisionTemp);
%     eddySizeDivision=cellfun(@(x) mean(x),eddySizeDivisionTemp);
%     eddyOWDivision=cellfun(@(x) mean(x),eddyOWDivisionTemp);
    

    % Visualization of dataset below
    figure,  
    ax4 = axes;

    eddySizeDivision=imgaussfilt(eddySizeDivision,1);
    eddySizeDivision=eddySizeDivision./length(FrameNum);
    h1 = uimage(ax4,xSize,ySize,eddySizeDivision',"CDataMapping","scaled");
    hold on 
    for i = 1:1:size(B,1)
        boundary=B{i};
        plot(x_val(boundary(:,1)), y_val(boundary(:,2)),'k','LineWidth',2);
        hold on
    end
    cb = colorbar();
    cb.Label.String='Average Eddy Size';
%     set(h1, "alphadata",~isnan(eddySizeDivision)');
    pbaspect([1,1,1]);
    title("Avereage Eddy Size for whole dataset("+int2str(patchSize*7)+"km x "+int2str(patchSize*7)+"km per patch)");
    ylabel('Lattude');
    xlabel('Longitude');
%     caxis([4,10]);
    set(ax4,'YDir','normal');
    

    
    figure,  
    ax5 = axes;
    eddyNumberDivision=imgaussfilt(eddyNumberDivision,1);
    eddyNumberDivision=eddyNumberDivision./length(FrameNum);
    h2 = uimage(ax5,xSize,ySize,eddyNumberDivision',"CDataMapping","scaled");
    hold on 
    for i = 1:1:size(B,1)
        boundary=B{i};
        plot(x_val(boundary(:,1)), y_val(boundary(:,2)),'k','LineWidth',2);
        hold on
    end
    cb = colorbar();
    cb.Label.String='Average Eddy Num';
%     set(h2, "alphadata",eddyNumberDivision');
    pbaspect([1,1,1]);
    title("Avereage Eddy Number for whole dataset("+int2str(patchSize)+"pixel x "+int2str(patchSize)+"pixel per patch)");
    ylabel('Lattude');
    xlabel('Longitude');
%     caxis([0,5]);
    set(ax5,'YDir','normal');
    
    figure,  
    ax6 = axes;   
    eddyOWDivision=imgaussfilt(eddyOWDivision,1);
    eddyOWDivision=eddyOWDivision./length(FrameNum);
    
    h3 = uimage(ax6,xSize,ySize,eddyOWDivision',"CDataMapping","scaled");
    hold on 
    for i = 1:1:size(B,1)
        boundary=B{i};
        plot(x_val(boundary(:,1)), y_val(boundary(:,2)),'k','LineWidth',2);
        hold on
    end
    cb = colorbar();
    cb.Label.String='Average Eddy Vorticity';
%     set(h3, "alphadata",~isnan(eddyOWDivision)');
    pbaspect([1,1,1]);
    title("Avereage Eddy Vorticity for whole dataset("+int2str(patchSize)+"pixel x "+int2str(patchSize)+"pixel per patch)");
    ylabel('Lattude');
    xlabel('Longitude');
%     caxis([0,0.006]);
    set(ax6,'YDir','normal');
    
        figure,  
    ax7 = axes;   
    eddyNumTimeSizeDivision=imgaussfilt(eddyNumTimeSizeDivision,1);
    eddyNumTimeSizeDivision=eddyNumTimeSizeDivision./length(FrameNum);
    
    h3 = uimage(ax7,xSize,ySize,eddyNumTimeSizeDivision',"CDataMapping","scaled");
    hold on 
    for i = 1:1:size(B,1)
        boundary=B{i};
        plot(x_val(boundary(:,1)), y_val(boundary(:,2)),'k','LineWidth',2);
        hold on
    end
    cb = colorbar();
    cb.Label.String='Average Eddy Size Times Nums';
%     set(h3, "alphadata",~isnan(eddyOWDivision)');
    pbaspect([1,1,1]);
    title("Avereage Eddy Size Times Nums for whole dataset("+int2str(patchSize)+"pixel x "+int2str(patchSize)+"pixel per patch)");
    ylabel('Lattude');
    xlabel('Longitude');
%     caxis([0,0.006]);
    set(ax7,'YDir','normal');
end

