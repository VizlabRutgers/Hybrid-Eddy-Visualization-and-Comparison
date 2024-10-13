function [] = eddyVisStat(clockEddy,counterclockEddy,eddyPath,srcData)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

close all;
x_val = double(ncread(srcData, "xh"));
y_val = double(ncread(srcData, "yh"));
z_val = double(ncread(srcData, "z_l"));

FrameNum=1;
startLoc = [1,1,1,FrameNum];
count = [length(x_val),length(y_val),length(z_val),1];
stride = [1,1,1,1];

startLoc_2D = [1,1,1,FrameNum];
count_2D = [length(x_val),length(y_val),1,1];
stride_2D = [1,1,1,1];


u_val = ncread(srcData, "u", startLoc_2D, count_2D, stride_2D);
v_val = ncread(srcData, "v",startLoc_2D, count_2D, stride_2D);

figure
quiver(x_val, y_val, u_val', v_val');
hold on

eddy_x_clock = clockEddy(:,1);
eddy_y_clock = clockEddy(:,2);
eddy_x_counterclock = counterclockEddy(:,1);
eddy_y_counterclock = counterclockEddy(:,2);





clock_eddy_scatter = scatter(eddy_x_clock, eddy_y_clock, 8, "filled");
counterclock_eddy_scatter = scatter(eddy_x_counterclock, eddy_y_counterclock, 8, "filled");
daspect([1 1 1]);
legend([clock_eddy_scatter, counterclock_eddy_scatter], "clockwise", "counterclockwise");
xlabel("longitude");
ylabel("latitude");
title("Eddy spatial distribution");

radius_clock_maximum = max(clockEddy(:, 13));
radius_counterclock_maximum = max(counterclockEddy(:, 13));
radius_clock = hist(clockEddy(:, 13),radius_clock_maximum)./length(clockEddy(:, 13));
radius_counterclock = hist(counterclockEddy(:, 13),radius_counterclock_maximum)./length(counterclockEddy(:, 13));


figure,
clock_radius_hist = plot(1:1:radius_clock_maximum+2, [0,0,radius_clock]);
hold on
counterclock_radius_hist = plot(1:1:radius_counterclock_maximum+2,[0,0,radius_counterclock]);
legend([clock_radius_hist, counterclock_radius_hist], "clockwise", "counterclockwise");
xlabel("radius (pixel)")
ylabel("Distribution probability");
title("radius distribution");
end