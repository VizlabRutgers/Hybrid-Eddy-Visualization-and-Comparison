function [evaluationResult_all] = eddyInterpolation(pathIndex,G,eddyPath,eddyPathHistory,dataFilePath,srcData, property, eddyNumber)
%EDDYINTERPOLATION 
close all;

% fh1 = figure();
% 
% ax1=subplot(1,3,1);
% ax2 = subplot(1,3,2);
% ax3 = subplot(1,3,3);

evaluationResult_all = [];


x_val = double(ncread(srcData, property.x));
y_val = double(ncread(srcData, property.y));
z_val = double(ncread(srcData, property.z));



% FrameNum=1;
% startLoc = [1,1,1,FrameNum];
% count = [length(x_val),length(y_val),length(z_val),1];
% stride = [1,1,1,1];
% 
% startLoc_2D = [1,1,1,FrameNum];
% count_2D = [length(x_val),length(y_val),1,1];
% stride_2D = [1,1,1,1];
% 
% 
% u_val = ncread(srcData, property.u, startLoc_2D, count_2D, stride_2D);
% v_val = ncread(srcData, property.v,startLoc_2D, count_2D, stride_2D);

% For each eddy
% for i = pathIndex
i = eddyNumber;

eddyWholePath = eddyPathHistory{i};

wholeEddyDepth = {};
wholecenterX = {};
wholecenterY = {};

% Each frame in the eddy history
for eddyFrameIndex = 1:1:size(eddyWholePath,1)
    FrameNum = eddyWholePath(eddyFrameIndex, 15);

    if(eddyWholePath(eddyFrameIndex,14)==1)
        data = load(dataFilePath+"Seperated Structures/clockwise/Frame_"+num2str(FrameNum)+"_eddy_"+num2str(eddyWholePath(eddyFrameIndex,16))+"_statistic.uocd");
    elseif(eddyWholePath(eddyFrameIndex,14)==0)
        data = load(dataFilePath+"Seperated Structures/counterclockwise/Frame_"+num2str(FrameNum)+"_eddy_"+num2str(eddyWholePath(eddyFrameIndex,16))+"_statistic.uocd");
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

    depth = unique(data(:,5));
    [~, depthIndexTotal] = ismember(round(depth), round(z_val));

    minX=min(data(:,3));
    maxX=max(data(:,3));
    minZ=min(data(:,5));
    maxZ=max(data(:,5));

%         if(numel(depthIndexTotal) < 3)
% 
%             wholecenterX{1, FrameNum} = 0;
%             wholecenterY{1, FrameNum} = 0;
%             wholeEddyDepth{1, FrameNum} = 0;
%     
%             wholecenterX{2, FrameNum} = 0;
%             wholecenterY{2, FrameNum} = 0;
%             wholeEddyDepth{2, FrameNum} = 0;
% 
%             continue;
%         end
    
    % Get the boundary of this eddy structure
%         bound = boundary(x,y,z,0.9);
%         localOb = trisurf(bound,x,y,z,sqrt(u_obj.^2+v_obj.^2),'EdgeColor','none', 'FaceAlpha', '0.9','Parent', ax1);
%         localOb.SpecularExponent = 200;
%         localOb.AmbientStrength = 0.8;
%         hold(ax1,'off');
%         set(ax1, "ZDir", "reverse");
%         xlim1 = xlim(ax1);
%         ylim1 = ylim(ax1);


    shape_left=[];
    shape_right=[];
    radiusRecordTemp=[];
    centerXRecord = [];
    centerYRecord = [];

    % Get the center & depth vs. radius of eddy in one frame
    for depthIndex = 1:1:length(depthIndexTotal)
    
        dataInDepth = data(data(:,5)==depth(depthIndex),:);

        centerX = dataInDepth(1,1);
        centerY = dataInDepth(1,2);
        x = dataInDepth(:,3);
        y = dataInDepth(:,4);
        zInDepth = dataInDepth(:,5);
        ow = dataInDepth(:,6);
        u_obj = dataInDepth(:,7);
        v_obj = dataInDepth(:,8);
        temp = dataInDepth(:,10);
        salt = dataInDepth(:,11);

        radius = (max(x) - min(x))/2;

        shape_left = [shape_left; min(x),depth(depthIndex)];
        shape_right = [shape_right; max(x),depth(depthIndex)];
        radiusRecordTemp = [radiusRecordTemp; radius];
        centerXRecord = [centerXRecord, centerX];
        centerYRecord = [centerYRecord, centerY];
    end
    
%         plot(ax2,radiusRecordTemp,'-o');
% 
%         x2 = xlabel(ax2,'depth');
%         y2 = ylabel(ax2,'radius (distance in degree)');
% 
%         plot3(ax3,centerXRecord, centerYRecord, unique(z),'-o');
%         x1 = xlabel(ax1,'Longitude');
%         y1 = ylabel(ax1,'Latitude');
%         
%         grid(ax3, "on");
%         set(ax3, "ZDir", "reverse");
%         x3 = xlabel(ax3,'Longitude');
%         y3 = ylabel(ax3,'Latitude');
%         z3 = zlabel(ax3,'depth');
%         xlim([xlim1(1) xlim1(2)])
%         ylim([ylim1(1) ylim1(2)])

    wholecenterX{1, FrameNum} = centerXRecord;
    wholecenterY{1, FrameNum} = centerYRecord;
    wholeEddyDepth{1, FrameNum} = radiusRecordTemp;

    wholecenterX{2, FrameNum} = numel(depthIndexTotal);
    wholecenterY{2, FrameNum} = numel(depthIndexTotal);
    wholeEddyDepth{2, FrameNum} = numel(depthIndexTotal);
end




% DTWDTWDTWDTWDTW
% DTWDTWDTWDTWDTW
% DTWDTWDTWDTWDTW
% DTWDTWDTWDTWDTW
% DTWDTWDTWDTWDTW




for eddyFrameIndex = 1:1:(size(eddyWholePath,1)-2)
    startFrame = eddyWholePath(eddyFrameIndex, 15);   



    maxStride = 6;

    if((startFrame + maxStride) > eddyWholePath(end, 15))
        maxStride = eddyWholePath(end, 15) - startFrame;
    end


%     startFrame = 19;
%     originData_stride = 2;

    for originData_stride = 2:1:maxStride

%         if(startFrame == 57 && originData_stride == 2)
%             test=0;
%         end

        evaluationResult = interpolation_internal_alter(wholecenterX, wholecenterY, wholeEddyDepth, startFrame, originData_stride,z_val);

        cell_eddyNumber = cell(size(evaluationResult,1),1);
        cell_eddyNumber(:) = {eddyNumber};

        cell_interpolationStride = cell(size(evaluationResult,1),1);
        cell_interpolationStride(:) = {startFrame+originData_stride};



        evaluationResult = [cell_eddyNumber,evaluationResult(:,6),evaluationResult(:,5),evaluationResult(:,7),cell_interpolationStride,evaluationResult(:,8),evaluationResult(:,1:4),evaluationResult(:,13),evaluationResult(:,9:12),evaluationResult(:,14:17),evaluationResult(:,18)];
        

        
        evaluationResult_all = [evaluationResult_all;evaluationResult];
    end
    


end



%% physical attributes interpolation

% physicalAttrInterpolation();

end



