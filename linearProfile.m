function [outputArg1,outputArg2] = linearProfile(srcData, property, ...
    x_lower, x_upper, y_lower, y_upper, FrameNum, allEddy_hybrid,...
    eddySummary_oldWindingAngle, eddySummary_newWindingAngle, plotHybridOnly_flag, fixed_linearProfile_flag)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
close all;

%% preprocessing
x_val = double(ncread(srcData, property.x));
y_val = double(ncread(srcData, property.y));
z_val = double(ncread(srcData, property.z));

startLoc = [1,1,1,FrameNum];
count = [length(x_val),length(y_val),length(z_val),1];
stride = [1,1,1,1];

startLoc_surface = [1,1,1,FrameNum];
count_surface = [length(x_val),length(y_val),1,1];
stride_surface = [1,1,1,1];

startLoc_2D = [1,1,FrameNum];
count_2D = [length(x_val),length(y_val),1];
stride_2D = [1,1,1];

u_val = ncread(srcData, property.u, startLoc_surface, count_surface, stride_surface);
v_val = ncread(srcData, property.v,startLoc_surface, count_surface, stride_surface);
%     w_val = ncread(srcData, "W",startLoc, count, stride);
eta_val = ncread(srcData,property.eta,startLoc_2D, count_2D, stride_2D);
velocity_mag_val = (sqrt(u_val.^2 +v_val.^2));

xSize = length(x_val);
ySize = length(y_val);
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
% figure,imagesc(x_val,y_val,omega');
% colorbar();
% caxis([0,0.02])
% set(gca, "YDir", "normal");

% figure,imagesc(x_val,y_val,sqrt(u_val.^2+v_val.^2)');
% colorbar();
% caxis([0,0.4])
% set(gca, "YDir", "normal");
% ow_val = ow_val;

%% find the center for each approach

% old winding angle
eddyHistory_old_currentFrame = eddySummary_oldWindingAngle{FrameNum}{1};
currentCenter_X_old = eddyHistory_old_currentFrame{1};
currentCenter_Y_old = eddyHistory_old_currentFrame{2};
currentoffset_X = eddyHistory_old_currentFrame{4};
currentoffset_Y = eddyHistory_old_currentFrame{5};

[~,oldIndex] = min((currentCenter_X_old-(x_upper+x_lower)/2).^2 +...
    (currentCenter_Y_old-(y_upper+y_lower)/2).^2);

% new winding angle
eddyHistory_new_currentFrame = eddySummary_newWindingAngle{FrameNum}{1};
currentCenter_X_new = eddyHistory_new_currentFrame{1};
currentCenter_Y_new = eddyHistory_new_currentFrame{2};
currentboundary_X = eddyHistory_new_currentFrame{6};
currentboundary_Y = eddyHistory_new_currentFrame{7};

[~,newIndex] = min((currentCenter_X_new-(x_upper+x_lower)/2).^2 +...
    (currentCenter_Y_new-(y_upper+y_lower)/2).^2);

% hybrid
currEddydata = allEddy_hybrid(allEddy_hybrid(:,15)==FrameNum,:);

[~,hybridIndex] = min((currEddydata(:,1)-(x_upper+x_lower)/2).^2 +...
    (currEddydata(:,2)-(y_upper+y_lower)/2).^2);

hybrid_centerX = currEddydata(hybridIndex,1);
hybrid_centerY = currEddydata(hybridIndex,2);
radius = currEddydata(hybridIndex,13);

%% Calculate profile line
[~,x_lowerIndex] = min(abs(x_lower-x_val));
[~,x_upperIndex] = min(abs(x_upper-x_val));
[~,y_lowerIndex] = min(abs(y_lower-y_val));
[~,y_upperIndex] = min(abs(y_upper-y_val));

if(fixed_linearProfile_flag == 1)
    % Do nothing
    point1(1) = x_lower;
    point1(2) = hybrid_centerY;
    point2(1) = x_upper;
    point2(2) = hybrid_centerY;

elseif(fixed_linearProfile_flag == 0)
    % Switch the point to the center of new & old winding angle centers
    point1(1) = currentCenter_X_old(oldIndex);
    point1(2) = currentCenter_Y_old(oldIndex);
    point2(1) = currentCenter_X_new(newIndex);
    point2(2) = currentCenter_Y_new(newIndex);
end

[~,pt1_x] = min(abs(point1(1)-x_val));
[~,pt2_x] = min(abs(point2(1)-x_val));
[~,pt1_y] = min(abs(point1(2)-y_val));
[~,pt2_y] = min(abs(point2(2)-y_val));

pt1_x = pt1_x-x_lowerIndex;
pt2_x = pt2_x-x_lowerIndex;
pt1_y = pt1_y-y_lowerIndex;
pt2_y = pt2_y-y_lowerIndex;


M = x_upperIndex- x_lowerIndex+1;
N = y_upperIndex- y_lowerIndex+1;

% Coord range (x,y) of the line
scale = max(M,N)*10; 
x_line = [pt1_x pt2_x] + [-scale scale]*(pt2_x-pt1_x);
y_line = [pt1_y pt2_y] + [-scale scale]*(pt2_y-pt1_y);

xv = [1 M M 1 1];
yv = [1 1 N N 1];

[xi, yi] = polyxpoly(x_line, y_line, xv, yv, 'unique');
dx = xi(1)-xi(2); 
dy = yi(1)-yi(2);
L = round(sqrt(dx^2 + dy^2));

figure_overall=figure();


ax_overall = axes(figure_overall);
imagesc(ax_overall,x_val, y_val,sqrt(u_val.^2 + v_val.^2)');
cb = colorbar();
set(ax_overall, "YDir", "normal");
xlabel("longitude");
ylabel("latitude");
cb.Label.String = "Velocity magnitude";
daspect([1 1 1])
ax_overall.FontSize = 24;


%% ----------Old Winding Angle Method---------




oldWindingAngleX =currentCenter_X_old(oldIndex)+currentoffset_X{oldIndex};
oldWindingAngleY =currentCenter_Y_old(oldIndex)+currentoffset_Y{oldIndex};

[~,oldWindingAngleX] = min(abs(oldWindingAngleX-x_val));
[~,oldWindingAngleY] = min(abs(oldWindingAngleY-y_val));

oldWindingAngleX = oldWindingAngleX-x_lowerIndex;
oldWindingAngleY = oldWindingAngleY-y_lowerIndex;

[x_oldWinding, y_oldWinding] = polyxpoly(x_line, y_line, oldWindingAngleX, oldWindingAngleY, 'unique');

for k=1:numel(x_oldWinding)
    t = sqrt(((x_oldWinding(k)-xi(2)+1)^2 + (y_oldWinding(k)-yi(2)+1)^2)/(dx^2+dy^2));
    boundaryOnProfile_oldWindingAngle(k) = t*L;
end


[~, currentCenter_X_old_single] = min(abs(currentCenter_X_old(oldIndex)-x_val));
[~, currentCenter_Y_old_single] = min(abs(currentCenter_Y_old(oldIndex)-y_val));

t = sqrt(((currentCenter_X_old_single-x_lowerIndex-xi(2)+2)^2 + (currentCenter_Y_old_single-y_lowerIndex-yi(2)+2)^2)/(dx^2+dy^2));
centerOnProfile_oldWindingAngle = t*L;


    %% ----------New Winding Angle Method---------




[~,newWindingAngleX] = min(abs(currentboundary_X{newIndex}'-x_val));
[~,newWindingAngleY] = min(abs(currentboundary_Y{newIndex}'-y_val));

newWindingAngleX = newWindingAngleX-x_lowerIndex;
newWindingAngleY = newWindingAngleY-y_lowerIndex;

[x_newWinding, y_newWinding] = polyxpoly(x_line, y_line, newWindingAngleX, newWindingAngleY, 'unique');

for k=1:numel(x_newWinding)
    t = sqrt(((x_newWinding(k)-xi(2)+1)^2 + (y_newWinding(k)-yi(2)+1)^2)/(dx^2+dy^2));
    boundaryOnProfile_newWindingAngle(k) = t*L;
end 

[~, currentCenter_X_new_single] = min(abs(currentCenter_X_new(newIndex)-x_val));
[~, currentCenter_Y_new_single] = min(abs(currentCenter_Y_new(newIndex)-y_val));

t = sqrt(((currentCenter_X_new_single-x_lowerIndex-xi(2)+2)^2 + (currentCenter_Y_new_single-y_lowerIndex-yi(2)+2)^2)/(dx^2+dy^2));
centerOnProfile_newWindingAngle = t*L;


%% ----------Hybrid Method---------


r = radius*0.04;  
pos = [hybrid_centerX,hybrid_centerY]; 
t=0:0.001:(2*pi);  
t=[t,0];
hybrid_boundary_X = pos(1)+r*sin(t);
hybrid_boundary_Y = pos(2)+r*cos(t);


[~,hybridBoundaryX] = min(abs(hybrid_boundary_X-x_val));
[~,hybridBoundaryY] = min(abs(hybrid_boundary_Y-y_val));

hybridBoundaryX = hybridBoundaryX-x_lowerIndex;
hybridBoundaryY = hybridBoundaryY-y_lowerIndex;

[x_newHybrid, y_newHybrid] = polyxpoly(x_line, y_line, hybridBoundaryX, hybridBoundaryY, 'unique');

for k=1:numel(x_newHybrid)
    t = sqrt(((x_newHybrid(k)-xi(2)+1)^2 + (y_newHybrid(k)-yi(2)+1)^2)/(dx^2+dy^2));
    boundaryOnProfile_hybrid(k) = t*L;
end 

[~, hybrid_centerX_single] = min(abs(hybrid_centerX-x_val));
[~, hybrid_centerY_single] = min(abs(hybrid_centerY-y_val));

t = sqrt(((hybrid_centerX_single-x_lowerIndex-xi(2)+2)^2 + (hybrid_centerY_single-y_lowerIndex-yi(2)+2)^2)/(dx^2+dy^2));
centerOnProfile_hybrid = t*L;
 

%% horizontal plane


ow_patch = ow_val(x_lowerIndex:x_upperIndex, y_lowerIndex:y_upperIndex);
figure, imagesc(x_val(x_lowerIndex:x_upperIndex), y_val(y_lowerIndex:y_upperIndex),...
    ow_patch');
% title("OW");
hold on
plot(x_val(round(xi+x_lowerIndex-1)), y_val(round(yi+y_lowerIndex)), "w--", "linewidth", 2.5);
hold on

if(plotHybridOnly_flag == 0)
    newWindingAngleCenter = plot(currentCenter_X_new(newIndex), currentCenter_Y_new(newIndex), "r.", 'MarkerSize', 40);
    hold on
    oldWindingAngleCenter = plot(currentCenter_X_old(newIndex), currentCenter_Y_old(newIndex), "k.", 'MarkerSize', 18);
    hold on
end
hybridCenter = plot(hybrid_centerX, hybrid_centerY, '.','color',[0.4940, 0.1840, 0.5560],...
    'MarkerSize', 28);
hold on
colormap(pink)
cb = colorbar();
cb.Label.String = "OW";
daspect([1 1 1]);
caxis([-0.01, 0.008]);
set(gca, "YDir", "normal");
hold on
hybrid_boundary= plot(hybrid_boundary_X,hybrid_boundary_Y,'-','color',[0.4940, 0.1840, 0.5560],...
    'linewidth',4);
if(plotHybridOnly_flag == 0)
    newWinding_boundary = plot(currentboundary_X{newIndex}, currentboundary_Y{newIndex}, ...
        'r','linewidth',4);
    oldWinding_boundary = plot(currentCenter_X_old(oldIndex)+currentoffset_X{oldIndex},...
        currentCenter_Y_old(oldIndex)+currentoffset_Y{oldIndex}, 'k','linewidth',4);
%     lg = legend([newWinding_boundary,oldWinding_boundary,hybrid_boundary,newWindingAngleCenter,...
%         oldWindingAngleCenter,hybridCenter],...
%         "New winding angle boundary", "Old winding angle boundary","Hybrid boundary", ...
%         "New winding angle cener", "Old winding angle center",...
%         "Hybrid center");
else
%     lg = legend([hybrid_boundary, hybridCenter],...
%         "Hybrid boundary","Hybrid center");
end
% xlabel("longitude");
% ylabel("latitude");
% lg.FontSize = 10;
cb.FontSize = 18;
axis off
hold off

velmag_patch = velocity_mag_val(x_lowerIndex:x_upperIndex, y_lowerIndex:y_upperIndex);
figure, imagesc(x_val(x_lowerIndex:x_upperIndex), y_val(y_lowerIndex:y_upperIndex),...
    velmag_patch');
% title("velocity magnitude");
hold on
plot(x_val(round(xi+x_lowerIndex)-1), y_val(round(yi+y_lowerIndex)), "w--", "linewidth", 2.5);

if(plotHybridOnly_flag == 0)
    newWindingAngleCenter = plot(currentCenter_X_new(newIndex), currentCenter_Y_new(newIndex), "r.", 'MarkerSize', 40);
    hold on
    oldWindingAngleCenter = plot(currentCenter_X_old(newIndex), currentCenter_Y_old(newIndex), "k.", 'MarkerSize', 18);
    hold on
end
hold on
hybridCenter = plot(hybrid_centerX, hybrid_centerY, '.','color',[0.4940, 0.1840, 0.5560],...
    'MarkerSize', 28);
hold on
cb = colorbar();
cb.Label.String = "Velocity Magnitude";
daspect([1 1 1]);
caxis([0, 0.5]);
set(gca, "YDir", "normal");
hybrid_boundary= plot(hybrid_boundary_X,hybrid_boundary_Y,'-','color',[0.4940, 0.1840, 0.5560],...
    'linewidth',4);
if(plotHybridOnly_flag == 0)

    oldWinding_boundary = plot(currentCenter_X_old(oldIndex)+currentoffset_X{oldIndex},...
        currentCenter_Y_old(oldIndex)+currentoffset_Y{oldIndex}, 'k','linewidth',4);
    newWinding_boundary = plot(currentboundary_X{newIndex}, currentboundary_Y{newIndex}, ...
        'r','linewidth',4);
%     lg = legend([newWinding_boundary,oldWinding_boundary,hybrid_boundary,newWindingAngleCenter,...
%         oldWindingAngleCenter,hybridCenter],...
%         "New winding angle boundary", "Old winding angle boundary","Hybrid boundary", ...
%         "New winding angle cener", "Old winding angle center",...
%         "Hybrid center");
else
%     lg = legend([hybrid_boundary, hybridCenter],...
%         "Hybrid boundary","Hybrid center");
end
% xlabel("longitude");
% ylabel("latitude");
cb.FontSize = 18;
axis off
hold off

ssh_patch = eta_val(x_lowerIndex:x_upperIndex, y_lowerIndex:y_upperIndex);
figure, imagesc(x_val(x_lowerIndex:x_upperIndex), y_val(y_lowerIndex:y_upperIndex),...
    ssh_patch');
% title("SSH");
hold on
plot(x_val(round(xi+x_lowerIndex)-1), y_val(round(yi+y_lowerIndex)), "w--", "linewidth", 2.5);

if(plotHybridOnly_flag == 0)

    oldWindingAngleCenter = plot(currentCenter_X_old(newIndex), currentCenter_Y_old(newIndex), "k.", 'MarkerSize', 18);
    hold on
    newWindingAngleCenter = plot(currentCenter_X_new(newIndex), currentCenter_Y_new(newIndex), "r.", 'MarkerSize', 40);
    hold on
end
hybridCenter = plot(hybrid_centerX, hybrid_centerY, '.','color',[0.4940, 0.1840, 0.5560],...
    'MarkerSize', 28);
hold on
colormap(jet)
cb = colorbar();
cb.Label.String = "SSH";
daspect([1 1 1]);
set(gca, "YDir", "normal");
hybrid_boundary= plot(hybrid_boundary_X,hybrid_boundary_Y,'-','color',[0.4940, 0.1840, 0.5560],...
    'linewidth',4);
if(plotHybridOnly_flag == 0)

    oldWinding_boundary = plot(currentCenter_X_old(oldIndex)+currentoffset_X{oldIndex},...
        currentCenter_Y_old(oldIndex)+currentoffset_Y{oldIndex}, 'k','linewidth',4);
        newWinding_boundary = plot(currentboundary_X{newIndex}, currentboundary_Y{newIndex}, ...
        'r','linewidth',4);
%     lg = legend([newWinding_boundary,oldWinding_boundary,hybrid_boundary,newWindingAngleCenter,...
%         oldWindingAngleCenter,hybridCenter],...
%         "New winding angle boundary", "Old winding angle boundary","Hybrid boundary", ...
%         "New winding angle cener", "Old winding angle center",...
%         "Hybrid center");
else
%     lg = legend([hybrid_boundary, hybridCenter],...
%         "Hybrid boundary","Hybrid center");
end
% xlabel("longitude");
% ylabel("latitude");
cb.FontSize = 18;
axis off
hold off

%% linear profile

ow_profile_figure = figure();
ow_profile_axes = axes(ow_profile_figure);
ow_profile = improfile(ow_patch', [xi(2),xi(1)],[yi(2),yi(1)],L+1);
plot(ow_profile_axes,ow_profile,'linewidth',4);
title(ow_profile_axes,"OW profile");
% xlabel(ow_profile_axes,"distance on the profile");
ylabel(ow_profile_axes,"OW");
ylim([-0.02,0.02]);
ow_profile_axes.FontSize = 18;

if(plotHybridOnly_flag == 0)
    hold on
    xline(boundaryOnProfile_oldWindingAngle(1), 'k--', 'linewidth',4);
    hold on
    oldWindingAngleBoundary = xline(boundaryOnProfile_oldWindingAngle(2), 'k--', 'linewidth',4);
    xmarks = centerOnProfile_oldWindingAngle;
    ymarks = interp1(1:numel(ow_profile), ow_profile, xmarks, 'linearfilled');
    oldWindingAngleCenter = plot(xmarks, ymarks, 'k.', 'MarkerSize', 28);
    xline(boundaryOnProfile_newWindingAngle(1), 'r--', 'linewidth',4);
    hold on
    newWindingAngleBoundary = xline(boundaryOnProfile_newWindingAngle(2), 'r--', 'linewidth',4);
    hold on
    xmarks = centerOnProfile_newWindingAngle;  
    ymarks = interp1(1:numel(ow_profile), ow_profile, xmarks, 'linear');
    newWindingAngleCenter = plot(xmarks, ymarks, 'r.', 'MarkerSize', 40);
    hold on
end
xline(boundaryOnProfile_hybrid(1), '--', 'color',[0.4940, 0.1840, 0.5560],'linewidth',4);
hold on
hybridBoundary = xline(boundaryOnProfile_hybrid(2), '--', 'color',[0.4940, 0.1840, 0.5560],'linewidth',4);
xmarks = centerOnProfile_hybrid;  
ymarks = interp1(1:numel(ow_profile), ow_profile, xmarks, 'linear');
hybridCenter = plot(xmarks, ymarks, '.', 'color',[0.4940, 0.1840, 0.5560],'MarkerSize', 28);


if(plotHybridOnly_flag == 0)
%     lg = legend([newWindingAngleBoundary,oldWindingAngleBoundary,hybridBoundary,newWindingAngleCenter,oldWindingAngleCenter,...
%         hybridCenter],...
%         "New winding angle boundary", "Old winding angle boundary", "Hybrid boundary",...
%         "New winding angle cener", "Old winding angle center", "Hybrid center");
elseif(plotHybridOnly_flag == 1)
%     lg = legend([hybridBoundary,hybridCenter],...
%         "Hybrid boundary",...
%         "Hybrid center");
end
velmag_profile_figure = figure();
velmag_profile_axes = axes(velmag_profile_figure);
velmag_profile = improfile(velmag_patch', [xi(2),xi(1)],[yi(2),yi(1)],L+1);
plot(velmag_profile_axes,velmag_profile,'linewidth',4);
title(velmag_profile_axes,"Velocity magnitude profile");
% xlabel(velmag_profile_axes,"distance on the profile");
ylabel(velmag_profile_axes,"Velocity magnitude");
hold on
if(plotHybridOnly_flag == 0)
    hold on
    xline(boundaryOnProfile_oldWindingAngle(1), 'k--', 'linewidth',4);
    hold on
    oldWindingAngleBoundary = xline(boundaryOnProfile_oldWindingAngle(2), 'k--', 'linewidth',4);
    xmarks = centerOnProfile_oldWindingAngle;
    ymarks = interp1(1:numel(velmag_profile), velmag_profile, xmarks, 'linearfilled');
    oldWindingAngleCenter = plot(xmarks, ymarks, 'k.', 'MarkerSize', 28);
    xline(boundaryOnProfile_newWindingAngle(1), 'r--', 'linewidth',4);
    hold on
    newWindingAngleBoundary = xline(boundaryOnProfile_newWindingAngle(2), 'r--', 'linewidth',4);
    hold on
    xmarks = centerOnProfile_newWindingAngle;  
    ymarks = interp1(1:numel(velmag_profile), velmag_profile, xmarks, 'linear');
    newWindingAngleCenter = plot(xmarks, ymarks, 'r.', 'MarkerSize', 40);
    hold on
end
xline(boundaryOnProfile_hybrid(1), '--', 'color',[0.4940, 0.1840, 0.5560],'linewidth',4);
hold on
hybridBoundary = xline(boundaryOnProfile_hybrid(2), '--', 'color',[0.4940, 0.1840, 0.5560],'linewidth',4);
xmarks = centerOnProfile_hybrid;  
ymarks = interp1(1:numel(velmag_profile), velmag_profile, xmarks, 'linear');
hybridCenter = plot(xmarks, ymarks, '.', 'color',[0.4940, 0.1840, 0.5560],'MarkerSize', 28);
velmag_profile_axes.FontSize = 18;

if(plotHybridOnly_flag == 0)
%     lg = legend([newWindingAngleBoundary,oldWindingAngleBoundary,hybridBoundary,newWindingAngleCenter,oldWindingAngleCenter,...
%         hybridCenter],...
%         "New winding angle boundary", "Old winding angle boundary", "Hybrid boundary",...
%         "New winding angle cener", "Old winding angle center", "Hybrid center");
elseif(plotHybridOnly_flag == 1)
%     lg = legend([hybridBoundary,hybridCenter],...
%         "Hybrid boundary",...
%         "Hybrid center");
end
ssh_profile_figure = figure();
ssh_profile_axes = axes(ssh_profile_figure);
ssh_profile = improfile(ssh_patch', [xi(2),xi(1)],[yi(2),yi(1)],L+1);
plot(ssh_profile_axes,ssh_profile,'linewidth',4);
title(ssh_profile_axes,"SSH profile");
% xlabel(ssh_profile_axes,"distance on the profile");
ylabel(ssh_profile_axes,"SSH");
hold on
if(plotHybridOnly_flag == 0)
    hold on
    xline(boundaryOnProfile_oldWindingAngle(1), 'k--', 'linewidth',4);
    hold on
    oldWindingAngleBoundary = xline(boundaryOnProfile_oldWindingAngle(2), 'k--', 'linewidth',4);
    xmarks = centerOnProfile_oldWindingAngle;
    ymarks = interp1(1:numel(ssh_profile), ssh_profile, xmarks, 'linearfilled');
    oldWindingAngleCenter = plot(xmarks, ymarks, 'k.', 'MarkerSize', 28);
    xline(boundaryOnProfile_newWindingAngle(1), 'r--', 'linewidth',4);
    hold on
    newWindingAngleBoundary = xline(boundaryOnProfile_newWindingAngle(2), 'r--', 'linewidth',4);
    hold on
    xmarks = centerOnProfile_newWindingAngle;  
    ymarks = interp1(1:numel(ssh_profile), ssh_profile, xmarks, 'linear');
    newWindingAngleCenter = plot(xmarks, ymarks, 'r.', 'MarkerSize', 40);
    hold on
end
xline(boundaryOnProfile_hybrid(1), '--', 'color',[0.4940, 0.1840, 0.5560],'linewidth',4);
hold on
hybridBoundary = xline(boundaryOnProfile_hybrid(2), '--', 'color',[0.4940, 0.1840, 0.5560],'linewidth',4);
xmarks = centerOnProfile_hybrid;  
ymarks = interp1(1:numel(ssh_profile), ssh_profile, xmarks, 'linear');
hybridCenter = plot(xmarks, ymarks, '.', 'color',[0.4940, 0.1840, 0.5560],'MarkerSize', 28);
ssh_profile_axes.FontSize = 18;

if(plotHybridOnly_flag == 0)
%     lg = legend([newWindingAngleBoundary,oldWindingAngleBoundary,hybridBoundary,newWindingAngleCenter,oldWindingAngleCenter,...
%         hybridCenter],...
%         "New winding angle boundary", "Old winding angle boundary", "Hybrid boundary",...
%         "New winding angle cener", "Old winding angle center", "Hybrid center");
elseif(plotHybridOnly_flag == 1)
%     lg = legend([hybridBoundary,hybridCenter],...
%         "Hybrid boundary",...
%         "Hybrid center");
end
end