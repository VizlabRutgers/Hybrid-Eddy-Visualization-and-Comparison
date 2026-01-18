function [clkEddies,conclkEddies,EddyNumRecord] = eddyStatFunc(clockEddy,counterclockEddy,eddyHistory,srcData, property,durationLimits, frameLimits,z_val)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% Stat - part 1 - load data

% close all;
x_val = double(ncread(srcData, property.x));
y_val = double(ncread(srcData, property.y));
z_val = double(ncread(srcData, property.z));

FrameNum=1;
startLoc = [1,1,1,FrameNum];
count = [length(x_val),length(y_val),length(z_val),1];
stride = [1,1,1,1];

startLoc_2D = [1,1,1,FrameNum];
count_2D = [length(x_val),length(y_val),1,1];
stride_2D = [1,1,1,1];


u_val = ncread(srcData, property.u, startLoc_2D, count_2D, stride_2D);
v_val = ncread(srcData, property.v,startLoc_2D, count_2D, stride_2D);

distance_per_pixel = m_lldist([0,0], [y_val(1), y_val(2)]);







%% Stat - part 2 - filter

clkTempFrame=[];
conclkTempFrame=[];
clkTempmeanRadius=[];
conclkTempmeanRadius=[];
clkTempDistance=[];
conclkTempDistance=[];
alleddyTempFrame=[];
clkDepthTemp=[];
conclkDepthTemp=[];




clkeddyexistingFrame = [];
conclkeddyexistingFrame = [];
clkMovingDistance=[];
conclkMovingDistance=[];
clkmeanRadius=[];
conclkmeanRadius=[];
alleddyexistingFrame=[];
clkRadius=[];
conclkRadius=[];
clkDepth = [];
conclkDepth = [];




variable_names_types = [["id", "double"]; ...
			["startFrame", "double"]; ...
			["eddyexistingFrame", "double"]; ...
			["meanRadius", "double"]; ...
			["movingDistance", "double"];];
% Make table using fieldnames & value types from above
clkEddies = table('Size',[0,size(variable_names_types,1)],... 
	'VariableNames', variable_names_types(:,1),...
	'VariableTypes', variable_names_types(:,2));

conclkEddies = table('Size',[0,size(variable_names_types,1)],... 
	'VariableNames', variable_names_types(:,1),...
	'VariableTypes', variable_names_types(:,2));




NotRepeatFlag=0;

maxFrame = length(ncread(srcData, property.time));
clockwise_eddyMovingHistory= {};
counterclockwise_eddyMovingHistory= {};

EddyNumRecord = [];

for eddyIndex=1:1:length(eddyHistory)
    
%     if(eddyIndex == 17089) 
%         eddyIndex = 17089;
%     end

    eddyPathHistory=eddyHistory{eddyIndex};
%     thisEddytemp = arrayfun(@(x) eddyNodes(x), pathLine{1});
    
    thisEddy = eddyPathHistory(:,15:16)';
    
    if(size(thisEddy,2)<=durationLimits)
        continue;
    end

    if(thisEddy(1,1)<frameLimits.min)
        continue;
    end
    
    if(thisEddy(1,1)>frameLimits.max)
        continue;
    end




    % Seperate Eddy frame and index
%     thisEddytemp=out';
%     thisEddy=[];
%     thisEddy(1,1:length(thisEddytemp)/2)=thisEddytemp(2:2:end);
%     thisEddy(2,1:length(thisEddytemp)/2)=thisEddytemp(1:2:end);
    clockFlagLastTime=-1;
    
    thisEddyIndex = 1;
    
    previousCoord=[0,0];
    currentCoord=[0,0];
    innerCounter=0;

    clockwiseCheck=clockEddy((clockEddy(:,15)==thisEddy(1,1))&(clockEddy(:,16)==(thisEddy(2,1)-1)),:);
    conclockwiseCheck=counterclockEddy((counterclockEddy(:,15)==thisEddy(1,1))&(counterclockEddy(:,16)==(thisEddy(2,1)-1)),:);

%     eddyCheck=[];
%     for thisEddyIndex=1:1:size(thisEddy,2)  
%         % Eddy rotation check
%         NotRepeatFlag=0;
%         clocktempCheck=clockEddy((clockEddy(:,15)==thisEddy(1,thisEddyIndex))&(clockEddy(:,16)==(thisEddy(2,thisEddyIndex)-1)),:);
%         conclocktempCheck=counterclockEddy((counterclockEddy(:,15)==thisEddy(1,thisEddyIndex))&(counterclockEddy(:,16)==(thisEddy(2,thisEddyIndex)-1)),:);
%         if(~isempty(clocktempCheck))
%             eddyCheck(thisEddyIndex,:) = clocktempCheck;
%         elseif(~isempty(conclocktempCheck))
%             eddyCheck(thisEddyIndex,:) = conclocktempCheck;
%         end
%     end
    





%     if(mean(eddyCheck(:,1))<spaceLimit.x0) 
%         continue; 
%     end
%     if(mean(eddyCheck(:,1))>spaceLimit.x1) 
%         continue; 
%     end
%     if(mean(eddyCheck(:,2))<spaceLimit.y0) 
%         continue; 
%     end
%     if(mean(eddyCheck(:,2))>spaceLimit.y1) 
%         continue; 
%     end
%     if(mean(eddyCheck(:,13))<radiusLimit.lower)
%         continue;
%     end
% 
%     if(mean(eddyCheck(:,13))>radiusLimit.upper)
%         continue;
%     end
% 

%     if(~isempty(conclockwiseCheck)&&isempty(clockwiseCheck))      
%         if(mean(conclockwiseCheck(:,1))<spaceLimit.x0) 
%             continue; 
%         end
%         if(mean(conclockwiseCheck(:,1))>spaceLimit.x1) 
%             continue; 
%         end
%         if(mean(conclockwiseCheck(:,2))<spaceLimit.y0) 
%             continue; 
%         end
%         if(mean(conclockwiseCheck(:,2))>spaceLimit.y1) 
%             continue; 
%         end
%         if(mean(conclockwiseCheck(:,13))<radiusLimit.lower)
%             continue;
%         end
%     
%         if(mean(conclockwiseCheck(:,13))>radiusLimit.upper)
%             continue;
%         end
%     end

    clockwise_eddyMovingHistory_temp = [];
    counterclockwise_eddyMovingHistory_temp = [];
    thisEddyIndex = 1;

    while thisEddyIndex<=size(eddyPathHistory,1)  
        % Eddy rotation check
        innerCounter = innerCounter+1;
        NotRepeatFlag=0;
%         clockwiseCheck=clockEddy((clockEddy(:,15)==thisEddy(1,thisEddyIndex))&(clockEddy(:,16)==(thisEddy(2,thisEddyIndex)-1)),:);
%         conclockwiseCheck=counterclockEddy((counterclockEddy(:,15)==thisEddy(1,thisEddyIndex))&(counterclockEddy(:,16)==(thisEddy(2,thisEddyIndex))-1),:);
%         

        if(eddyPathHistory(thisEddyIndex,14)==1)
            clockwiseCheck = eddyPathHistory(thisEddyIndex,:);
        elseif(eddyPathHistory(thisEddyIndex,14)==0)
            conclockwiseCheck = eddyPathHistory(thisEddyIndex,:);
        else
            continue;
        end
        
        if(~isempty(clockwiseCheck)&&isempty(conclockwiseCheck))
            clockFlag=1;
        elseif(~isempty(conclockwiseCheck)&&isempty(clockwiseCheck))
            clockFlag=0;
        else
            thisEddyIndex=thisEddyIndex+1;
            continue;
        end
        
        
        if(clockFlag==clockFlagLastTime || clockFlagLastTime==-1)
            if(clockFlag==1)
                clkTempFrame=[clkTempFrame,1];
                [~,center_y] = ismembertol(clockwiseCheck(2),y_val,1e-4);
                clkTempmeanRadius=[clkTempmeanRadius,...
                    m_lldist([clockwiseCheck(1),clockwiseCheck(1)], ...
                    [y_val(center_y-clockwiseCheck(13)),...
                    y_val(center_y+clockwiseCheck(13))])/2];
                currentCoord=[clockwiseCheck(1),clockwiseCheck(2)];
                if(previousCoord(1)~=0)
                    clkTempDistance=[clkTempDistance,m_lldist([currentCoord(1),previousCoord(1)],[currentCoord(2),previousCoord(2)])];
                    conclkTempDistance=[conclkTempDistance,nan];

                end
                clockwise_eddyMovingHistory_temp = [clockwise_eddyMovingHistory_temp;currentCoord];
                clkDepthTemp=clockwiseCheck(5);

            else
                conclkTempFrame=[conclkTempFrame,1];
                [~,center_y] = ismembertol(conclockwiseCheck(2),y_val,1e-4);
                conclkTempmeanRadius=[conclkTempmeanRadius,...
                    m_lldist([conclockwiseCheck(1),conclockwiseCheck(1)], ...
                    [y_val(center_y-conclockwiseCheck(13)),...
                    y_val(center_y+conclockwiseCheck(13))])/2];
                currentCoord=[conclockwiseCheck(1),conclockwiseCheck(2)];
                if(previousCoord(1)~=0)
                    clkTempDistance=[clkTempDistance,nan];
                    conclkTempDistance=[conclkTempDistance,m_lldist([currentCoord(1),previousCoord(1)],[currentCoord(2),previousCoord(2)])];
                end
                counterclockwise_eddyMovingHistory_temp = [counterclockwise_eddyMovingHistory_temp;currentCoord];
                conclkDepthTemp=conclockwiseCheck(5);
            end
            previousCoord=currentCoord;
            alleddyTempFrame=[alleddyTempFrame,1];
            clockFlagLastTime = clockFlag;
        else
            
            if((length(clkTempFrame)>durationLimits)||(length(conclkTempFrame)>durationLimits))
            
                clkeddyexistingFrame=[clkeddyexistingFrame,length(clkTempFrame)];
                conclkeddyexistingFrame=[conclkeddyexistingFrame,length(conclkTempFrame)];
                clkmeanRadius=[clkmeanRadius,mean(clkTempmeanRadius)];
                conclkmeanRadius=[conclkmeanRadius,mean(conclkTempmeanRadius)];
                alleddyexistingFrame=[alleddyexistingFrame,length(alleddyTempFrame)];
                clkMovingDistance=[clkMovingDistance,sum(clkTempDistance)];
                conclkMovingDistance=[conclkMovingDistance,sum(conclkTempDistance)];
                clockwise_eddyMovingHistory= [clockwise_eddyMovingHistory;clockwise_eddyMovingHistory_temp];
                counterclockwise_eddyMovingHistory= [counterclockwise_eddyMovingHistory;counterclockwise_eddyMovingHistory_temp];
                clkRadius = [clkRadius, clkTempmeanRadius];
                conclkRadius = [conclkRadius, conclkTempmeanRadius];
                EddyNumRecord = [EddyNumRecord; eddyIndex];
                clkDepth = [clkDepth, mean(clkDepthTemp)];
                conclkDepth = [conclkDepth, mean(conclkDepthTemp)];
            end

            if(~isnan(mean(clkTempmeanRadius)))      
                clkEddyTemp = {eddyIndex, thisEddy(1,1),length(clkTempFrame),mean(clkTempmeanRadius),sum(clkTempDistance)};
                clkEddies = [clkEddies; clkEddyTemp];
            elseif(~isnan(mean(conclkTempmeanRadius)))
                conclkEddyTemp = {eddyIndex, thisEddy(1,1),length(conclkTempFrame),mean(conclkTempmeanRadius),sum(conclkTempDistance)};
                conclkEddies = [conclkEddies; conclkEddyTemp];
            end

            clkTempFrame=[];
            conclkTempFrame=[];
            clkTempmeanRadius=[];
            conclkTempmeanRadius=[];
            clkTempDistance=[];
            conclkTempDistance=[];
            alleddyTempFrame=[];
            clkDepthTemp=[];
            conclkDepthTemp=[];
            
            clockFlagLastTime = clockFlag;
            NotRepeatFlag=1;
            thisEddyIndex=thisEddyIndex+1;
            continue;
        end    
        thisEddyIndex=thisEddyIndex+1;
        
    end

    if(NotRepeatFlag==0)
        if((length(clkTempFrame)>=durationLimits)||(length(conclkTempFrame)>=durationLimits))
            clkeddyexistingFrame=[clkeddyexistingFrame,length(clkTempFrame)];
            conclkeddyexistingFrame=[conclkeddyexistingFrame,length(conclkTempFrame)];
            clkmeanRadius=[clkmeanRadius,mean(clkTempmeanRadius)];
            conclkmeanRadius=[conclkmeanRadius,mean(conclkTempmeanRadius)];
            alleddyexistingFrame=[alleddyexistingFrame,length(alleddyTempFrame)];
            clkMovingDistance=[clkMovingDistance,sum(clkTempDistance)];
            conclkMovingDistance=[conclkMovingDistance,sum(conclkTempDistance)];
            clockwise_eddyMovingHistory= [clockwise_eddyMovingHistory;clockwise_eddyMovingHistory_temp];
            counterclockwise_eddyMovingHistory= [counterclockwise_eddyMovingHistory;counterclockwise_eddyMovingHistory_temp];
            EddyNumRecord = [EddyNumRecord; eddyIndex];
            clkRadius = [clkRadius, clkTempmeanRadius];
            conclkRadius = [conclkRadius, conclkTempmeanRadius];
            clkDepth = [clkDepth, mean(clkDepthTemp)];
            conclkDepth = [conclkDepth, mean(conclkDepthTemp)];

        end

        if(~isnan(mean(clkTempmeanRadius)))      
            clkEddyTemp = {eddyIndex, thisEddy(1,1),length(clkTempFrame),mean(clkTempmeanRadius),sum(clkTempDistance)};
            clkEddies = [clkEddies; clkEddyTemp];
        elseif(~isnan(mean(conclkTempmeanRadius)))
            conclkEddyTemp = {eddyIndex, thisEddy(1,1),length(conclkTempFrame),mean(conclkTempmeanRadius),sum(conclkTempDistance)};
            conclkEddies = [conclkEddies; conclkEddyTemp];
        end

        clkTempFrame=[];
        conclkTempFrame=[];
        clkTempmeanRadius=[];
        conclkTempmeanRadius=[];
        clkTempDistance=[];
        conclkTempDistance=[];
        alleddyTempFrame=[];
        clkDepthTemp=[];
        conclkDepthTemp=[];
    end
    
    
    
end



clkeddyexistingFrame=nonzeros(clkeddyexistingFrame);
conclkeddyexistingFrame=nonzeros(conclkeddyexistingFrame);
clkmeanRadius=clkmeanRadius(~isnan(clkmeanRadius))';
conclkmeanRadius=conclkmeanRadius(~isnan(conclkmeanRadius))';
clkMovingDistance=clkMovingDistance(~isnan(clkMovingDistance));
conclkMovingDistance=conclkMovingDistance(~isnan(conclkMovingDistance));
alleddyTempFrame=nonzeros(alleddyTempFrame);
clkDepth = clkDepth(~isnan(clkDepth));
conclkDepth = conclkDepth(~isnan(conclkDepth));







%% Stat Visualization - part 1

clockEddy(clockEddy(:,15)<frameLimits.min,:)=[];
clockEddy(clockEddy(:,15)>frameLimits.max,:)=[];
counterclockEddy(counterclockEddy(:,15)<frameLimits.min,:)=[];
counterclockEddy(counterclockEddy(:,15)>frameLimits.max,:)=[];





figure,
subplot(1,2,1);
for i = 1:1:length(clockwise_eddyMovingHistory)
    clcMovingPath = clockwise_eddyMovingHistory{i,1};
    clcMovingPath = clcMovingPath - clcMovingPath(1,:);
    plot(clcMovingPath(:,1), clcMovingPath(:,2));
    hold on;
end
title('clockwise eddy moving history');
xlim([-5,5])
ylim([-5,5])
xlabel('longitude (\circ)');
ylabel('latitude (\circ)');






subplot(1,2,2);
for i = 1:1:length(counterclockwise_eddyMovingHistory)
    conclcMovingPath = counterclockwise_eddyMovingHistory{i,1};
    conclcMovingPath = conclcMovingPath - conclcMovingPath(1,:);
    plot(conclcMovingPath(:,1), conclcMovingPath(:,2));
    hold on;
end

title('counter-clockwise eddy moving history');
xlim([-5,5])
ylim([-5,5])
xlabel('longitude (\circ)');
ylabel('latitude (\circ)');


figure,
subplot(1,2,1);

clkMovingAngle = cellfun(@(x) atan2((x(end,2)-x(1,2)),(x(end,1)-x(1,1)+0.00000001)), clockwise_eddyMovingHistory, 'UniformOutput', true);
polarhistogram(clkMovingAngle,36)


title('clockwise eddy moving angle (start to end)');


subplot(1,2,2);
conclkMovingAngle = cellfun(@(x) atan2((x(end,2)-x(1,2)),(x(end,1)-x(1,1)+0.00000001)), counterclockwise_eddyMovingHistory, 'UniformOutput', true);
polarhistogram(conclkMovingAngle,36)

title('counter-clockwise eddy moving angle (start to end)');


figure
quiver(x_val, y_val, u_val', v_val');
hold on



% clockEddy(clockEddy(:,1)<spaceLimit.x0,:) = [];
% clockEddy(clockEddy(:,1)>spaceLimit.x1,:) = [];
% clockEddy(clockEddy(:,2)<spaceLimit.y0,:) = [];
% clockEddy(clockEddy(:,2)>spaceLimit.y1,:) = [];
% 
% counterclockEddy(counterclockEddy(:,1)<spaceLimit.x0,:) = [];
% counterclockEddy(counterclockEddy(:,1)>spaceLimit.x1,:) = [];
% counterclockEddy(counterclockEddy(:,2)<spaceLimit.y0,:) = [];
% counterclockEddy(counterclockEddy(:,2)>spaceLimit.y1,:) = [];

eddy_x_clock = clockEddy(:,1);
eddy_y_clock = clockEddy(:,2);
eddy_x_counterclock = counterclockEddy(:,1);
eddy_y_counterclock = counterclockEddy(:,2);




% clock_eddy_scatter = scatter(eddy_x_clock, eddy_y_clock, 8, "filled");
% counterclock_eddy_scatter = scatter(eddy_x_counterclock, eddy_y_counterclock, 8, "filled");
% daspect([1 1 1]);
% legend([clock_eddy_scatter, counterclock_eddy_scatter], "clockwise", "counterclockwise");
% xlabel("longitude");
% ylabel("latitude");
% title("Eddy spatial distribution");

subplot(1,2,1);
quiver(x_val, y_val, u_val', v_val');
hold on
clock_eddy_scatter = scatter(eddy_x_clock, eddy_y_clock, 8, "filled");
daspect([1 1 1]);
legend(clock_eddy_scatter, "clockwise");
xlabel("longitude");
ylabel("latitude");
title("Clockwise eddy spatial distribution");


subplot(1,2,2);
quiver(x_val, y_val, u_val', v_val');
hold on
counterclock_eddy_scatter = scatter(eddy_x_counterclock, eddy_y_counterclock, 8, "filled");
daspect([1 1 1]);
legend(counterclock_eddy_scatter,  "counterclockwise");
xlabel("longitude");
ylabel("latitude");
title("Counterclockwise eddy spatial distribution");



% for i=1:1:length(clockEddy)
%     [~,center_y(i)] = ismembertol(clockEddy(i,2),y_val,1e-4);
%     radius_clock(i) = m_lldist([clockEddy(i, 1),clockEddy(i, 1)], ...
%                     [y_val(center_y(i)-clockEddy(i, 13)),...
%                     y_val(center_y(i)+clockEddy(i, 13))])/2;
% end
% 
% for i=1:1:length(counterclockEddy)
%     [~,center_y(i)] = ismembertol(counterclockEddy(i,2),y_val,1e-4);
%     radius_counterclock(i) = m_lldist([counterclockEddy(i, 1),counterclockEddy(i, 1)], ...
%                     [y_val(center_y(i)-counterclockEddy(i, 13)),...
%                     y_val(center_y(i)+counterclockEddy(i, 13))])/2;
% end




[radius_clock_hist_result, radius_edges_clock]= histcounts(clkRadius,20);
[radius_counterclock_hist_result, radius_edges_counterclock] = histcounts(conclkRadius,20);


figure,
clock_radius_hist = plot([0,[radius_edges_clock(2:end)+radius_edges_clock(1:end-1)]/2], [0,radius_clock_hist_result]);
hold on
counterclock_radius_hist = plot([0,0,0,[radius_edges_counterclock(2:end)+radius_edges_counterclock(1:end-1)]/2],[0,0,0,radius_counterclock_hist_result]);
legend([clock_radius_hist, counterclock_radius_hist], "clockwise", "counterclockwise");
xlabel("radius (km)")
ylabel("Distribution probability");
title("radius distribution");

figure,
subplot(2,5,1);
[N,edges] = histcounts(clkeddyexistingFrame,'Normalization', 'probability');
histogram(clkeddyexistingFrame,'Normalization', 'probability','BinMethod','integers');
title('Clockwise eddy existing time');
xlabel('Living period (Frame)');
ylabel('probability');

subplot(2,5,6);
[N,edges] = histcounts(conclkeddyexistingFrame,'Normalization', 'probability');
histogram(conclkeddyexistingFrame,'Normalization', 'probability','BinMethod','integers');
title('Counterclockwise eddy existing time');
xlabel('Living period (Frame)');
ylabel('probability');

subplot(2,5,2);
[N,edges] = histcounts(clkMovingDistance,'Normalization', 'probability');
histogram(clkMovingDistance,40,'Normalization', 'probability');
title('Clockwise eddy moving distance');
xlabel('Distance (km)');
ylabel('probability');

subplot(2,5,7);
[N,edges] = histcounts(conclkMovingDistance,'Normalization', 'probability');
histogram(conclkMovingDistance,40,'Normalization', 'probability');
title('Counterclockwise eddy moving distance');
xlabel('Distance (km)');
ylabel('probability');

subplot(2,5,3);
scatter(clkmeanRadius, clkeddyexistingFrame, 'filled');
title('Average Radius vs existing time for each individual clockwise eddy over time');
xlabel('Radius (km)');
ylabel('Frames');

subplot(2,5,8);
scatter(conclkmeanRadius, conclkeddyexistingFrame, 'filled');
title('Average Radius vs existing time for each individual counter-clockwise eddy over time');
xlabel('Radius km)');
ylabel('Frames');

subplot(2,5,4);
scatter(clkmeanRadius, clkMovingDistance, 'filled');
title('Average Radius vs Moving distance for each individual clockwise eddy over time');
xlabel('Radius (km)');
ylabel('Distance (km)');

subplot(2,5,9);
scatter(conclkmeanRadius, conclkMovingDistance, 'filled');
title('Average Radius vs Moving distance for each individual counter-clockwise eddy over time');
xlabel('Radius (km)');
ylabel('Distance (km)');

subplot(2,5,5);
scatter(clkmeanRadius, clkDepth, 'filled');
title('Average Radius vs Moving distance for each individual clockwise eddy over time');
xlabel('Radius (km)');
ylabel('Depth (m)');

subplot(2,5,10);
scatter(conclkmeanRadius, conclkDepth, 'filled');
title('Average Radius vs Moving distance for each individual counter-clockwise eddy over time');
xlabel('Radius (km)');
ylabel('Depth (m)');


%% Stat Visualization - part 2
figure,
subplot(2,3,1);
boxchart(clkRadius);
hold on
plot(mean(clkRadius),'-o');
title("clockwise eddy radius distribution");

subplot(2,3,4);
boxchart(conclkRadius);
hold on
plot(mean(conclkRadius),'-o');
title("counterclockwise eddy radius distribution");

subplot(2,3,2);
boxchart(clkeddyexistingFrame);
hold on
plot(mean(clkeddyexistingFrame),'-o');
title("clockwise eddy exisinting frame distribution");

subplot(2,3,5);
boxchart(conclkeddyexistingFrame);
hold on
plot(mean(conclkeddyexistingFrame),'-o');
title("counterclockwise eddy exisinting frame distribution");

subplot(2,3,3);
boxchart(clkMovingDistance);
hold on
plot(mean(clkMovingDistance),'-o');
title("clockwise moving distance distribution");

subplot(2,3,6);
boxchart(conclkMovingDistance);
hold on
plot(mean(conclkMovingDistance),'-o');
title("counterclockwise moving distance distribution");


end

