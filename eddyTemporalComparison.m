function [] = eddyTemporalComparison(srcData, property, eddyIndex, eddyPathHistory,dataFilePath,stretchMode)
%UNTITLED 



% interpolationStartFrame = min(cell2mat(interpolationRecord(:,3)));
% interpolationEndFrame = max(cell2mat(interpolationRecord(:,5)));
% 
% interpolationRecord = interpolationRecord(cell2mat(interpolationRecord(:,2)) == stride4Test, :);
% 
% interpolationFrameIndex = interpolationStartFrame:stride4Test:interpolationEndFrame;
% interpolationRecord_filtered = interpolationRecord(ismember(cell2mat(interpolationRecord(:,3)),interpolationFrameIndex),:);


x_val = double(ncread(srcData, property.x));
y_val = double(ncread(srcData, property.y));
z_val = double(ncread(srcData, property.z));


hybridFlag = 0;
if(isfolder(dataFilePath+"Seperated Structures/clockwise/"))
    hybridFlag = 1;
end
    

fh2 = figure();
ax2=axes(fh2);


ax2.FontSize=24;

camlight(ax2);
lighting(ax2, 'flat');


hold(ax2,'on');
set(ax2,'ZDir','reverse');
view(3)



%% old winding angle

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
    

    % for index=1:1:size(eddyHistory,1)
    for index=1:1:size(eddyHistory,1)
        % Get each eddy in this eddy path history at this frame 
        % Plot the boundary of this eddy in 2D

        if(hybridFlag == 1)
            if(eddyHistory(index,14) == 1)
                data = load(dataFilePath+"Seperated Structures/clockwise/Frame_"+num2str(eddyHistory(index,15))+"_eddy_"+num2str(eddyHistory(index,16))+"_statistic.uocd");
            else
                data = load(dataFilePath+"Seperated Structures/counterclockwise/Frame_"+num2str(eddyHistory(index,15))+"_eddy_"+num2str(eddyHistory(index,16))+"_statistic.uocd");
            end
        else
            data = load(dataFilePath+"Seperated Structures/Frame_"+num2str(eddyHistory(index,15))+"_eddy_"+num2str(eddyHistory(index,16))+"_statistic.uocd");
        end
            
        horizontaloffset = 0;
        background_option = 0;
        eddyIndividual3DVis_inner(data,ax2,z_val, stretchMode, index, srcData, property,horizontaloffset, background_option);
        hold on
    end      
end

%% new widing angle
% for differentEddyIndex=eddyIndex_new
%     
% 
%     % Read eddy history
%     eddyHistory_new = eddyPathHistory_new{differentEddyIndex};
%     eddyHistory_new=sortrows(eddyHistory_new,15);
% 
% 
%     
%     % Extract the sub-graph of the eddy path
% 
%     
%     
%     
%     minX=inf;
%     maxX=0;
%     minY=inf;
%     maxY=0;
%     traceHistory=[];
%     
%     
%     
%     minX=min(eddyHistory_new(:,1));
%     maxX=max(eddyHistory_new(:,1));
%     minY=min(eddyHistory_new(:,2));
%     maxY=max(eddyHistory_new(:,2));
%     
%     
%     % for index=1:1:size(eddyHistory,1)
%     for index=1:1:10
%         % Get each eddy in this eddy path history at this frame 
%         % Plot the boundary of this eddy in 2D
% 
% 
%         data = load(dataFilePath_new+"Seperated Structures/Frame_"+num2str(eddyHistory_new(index,15))+"_eddy_"+num2str(eddyHistory_new(index,16))+"_statistic.uocd");
% 
%         horizontaloffset = 10;
%         eddyIndividual3DVis_inner(data,ax2,z_val, stretchMode, index, srcData, property, horizontaloffset);
%         hold on
%     end      
% end

set(ax2,'ZDir','reverse');
set(ax2, 'YDir', 'normal');
view(3)
daspect([1,1,100])    
ylim([-2,15]);
xlim([-4,26]);
zlim([0,1000]);
colormap(jet(5));
cb = colorbar();
cb.Label.String = "Velocity Magnitude";
caxis([0,0.5])
xlabel("Timestep");
ylabel("Actual size in Degree");
zlabel("Depth");
grid on
title("winding angle eddy detection results over 10 timesteps");

xticks(0:6:54);
xticklabels({'1','2','3','4','5','6','7','8','9','10'});

% grid on;
% set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
% cb.Label.String = "Velocity Magnitude";


% Remove y axis numbering




% lightangle(0,30);
lightangle(45,-45);
lighting(ax2, 'gouraud');
% lighting(ax2, 'flat');

% delete(findall(gcf, 'Type', 'light'));

end

