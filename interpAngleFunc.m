function [final_attr_interp] = interpAngleFunc(AttrInPrev, AttrInNext, AttrGroundtruth,originData_stride,interpolateIndex,interp_depth)
%INTERPFUNC Summary of this function goes here
%   Detailed explanation goes here
% Do the interpolation both frames to make sure they have
% same number of sample points (though they might have different depth)

deeperLayer = max(length(AttrInPrev), length(AttrInNext));
 

thetaPrev = deg2rad(AttrInPrev);
thetaNext = deg2rad(AttrInNext);


vPrev = [cos(thetaPrev), sin(thetaPrev)];
vNext = [cos(thetaNext), sin(thetaNext)]; 


xPrev_interp = interp1(1:length(AttrInPrev), vPrev(:,1)', 1:((length(AttrInPrev)-1)/(deeperLayer-1)):length(AttrInPrev))'; % 2 x N
yPrev_interp = interp1(1:length(AttrInPrev), vPrev(:,2)', 1:((length(AttrInPrev)-1)/(deeperLayer-1)):length(AttrInPrev))'; % 2 x N
xNext_interp = interp1(1:length(AttrInNext), vNext(:,1)', 1:(length(AttrInNext)-1)/(deeperLayer-1):length(AttrInNext))'; % 2 x N
yNext_interp = interp1(1:length(AttrInNext), vNext(:,2)', 1:(length(AttrInNext)-1)/(deeperLayer-1):length(AttrInNext))'; % 2 x N



w1 = (originData_stride - interpolateIndex) / originData_stride;
w2 = interpolateIndex / originData_stride;


x_interp = xPrev_interp * w1 + xNext_interp * w2;
y_interp = yPrev_interp * w1 + yNext_interp * w2;


norms = sqrt(x_interp.^2 + y_interp.^2);
x_interp = x_interp ./ norms;
y_interp = y_interp ./ norms;

interp_angles = rad2deg(atan2(y_interp, x_interp));
interp_angles(interp_angles < 0) = interp_angles(interp_angles < 0) + 360;

delta = mod((AttrGroundtruth(1) - interp_angles(1) + 180), 360) - 180;
correction = linspace(delta, 0, length(interp_angles));
interp_angles = mod(interp_angles + correction', 360);

final_attr_interp = interp1(interp_depth, interp_angles, 1:1:fix(interp_depth(end)));


    
end

