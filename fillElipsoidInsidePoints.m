function [x_inside,y_inside] = fillElipsoidInsidePoints(xc,yc,boundary_x, boundary_y, longAxis, shortAxis,theta,x_val,y_val)
%FILLELIPSOIDINSIDEPOINTS Summary of this function goes here
%   Detailed explanation goes here
xmin = min(boundary_x);
xmax = max(boundary_x);
ymin = min(boundary_y);
ymax = max(boundary_y);





x_val_minInd = find(x_val<xmin, 1, 'last');
x_val_maxInd = find(x_val>xmax, 1, 'first');
y_val_minInd = find(y_val<ymin, 1, 'last');
y_val_maxInd = find(y_val>ymax, 1, 'first');


[x_sub, y_sub] = meshgrid(x_val(x_val_minInd:x_val_maxInd), y_val(y_val_minInd:y_val_maxInd));






% 移动到椭圆中心坐标系
x_shift = x_sub - xc;
y_shift = y_sub - yc;

% 坐标旋转（逆旋转）
x_tmp = cosd(theta) * x_shift + sind(theta) * y_shift;
y_tmp = -sind(theta) * x_shift + cosd(theta) * y_shift;

% 使用标准椭圆公式判断
ellipse_mask = (x_tmp.^2 / longAxis^2 + y_tmp.^2 / shortAxis^2) <= 1;
x_inside = x_sub(ellipse_mask);
y_inside = y_sub(ellipse_mask);





end

