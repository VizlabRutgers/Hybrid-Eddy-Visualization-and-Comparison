function [outputArg1,outputArg2] = generateClip(srcData, property)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
x_val = double(ncread(srcData, property.x));
y_val = double(ncread(srcData, property.y));
z_val = double(ncread(srcData, property.z));

FrameNum=1;
startLoc = [338,1,1,FrameNum];
count = [401-338,76-1,length(z_val),1];
stride = [1,1,1,1];

startLoc_2D = [1,1,FrameNum];
count_2D = [length(x_val),length(y_val),1];
stride_2D = [1,1,1];


u_val = ncread(srcData, property.u, startLoc, count, stride);
v_val = ncread(srcData, property.v,startLoc, count, stride);
w_val = ncread(srcData, property.w,startLoc, count, stride);


end

