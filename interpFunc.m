function [final_attr_interp] = interpFunc(AttrInPrev, AttrInNext, AttrGroundtruth,originData_stride,interpolateIndex,interp_depth)
%INTERPFUNC Summary of this function goes here
%   Detailed explanation goes here
% Do the interpolation both frames to make sure they have
% same number of sample points (though they might have different depth)

deeperLayer = max(length(AttrInPrev), length(AttrInNext));
 
AttrInPrev_interp = interp1(1:1:length(AttrInPrev), AttrInPrev,1:((length(AttrInPrev)-1)/(deeperLayer-1)):length(AttrInPrev));
AttrInNext_interp = interp1(1:1:length(AttrInNext), AttrInNext, 1:(length(AttrInNext)-1)/(deeperLayer-1):length(AttrInNext));



% Associate the interpolated depth with corresponding depth,
% the deeper frame will keep unchanged
attr_depth_previousFrame = AttrInPrev_interp;
attr_depth_nextFrame = AttrInNext_interp;


% Do the interpolation based on the stride between the two
% frames
% Consider using the distance between two frames
attr_interpolationResult = attr_depth_previousFrame * (originData_stride-interpolateIndex) / originData_stride + attr_depth_nextFrame * interpolateIndex / originData_stride;

% adjust the interpolation result by the surface groundtruth
correction = linspace((AttrGroundtruth(1) - attr_interpolationResult(1)), 0, length(attr_interpolationResult));
attr_interpolationResult = attr_interpolationResult + correction;



final_attr_interp = interp1(interp_depth, attr_interpolationResult, 1:1:fix(interp_depth(end)));



    %             centerX_firstFrame_interp = interpolateData(centerX_firstFrame - centerX_firstFrame(1), interpolated_layer);
    %             centerY_firstFrame_interp = interpolateData(centerY_firstFrame - centerY_firstFrame(1), interpolated_layer);
    %     
    %             centerX_nextFrame_interp = interpolateData(centerX_nextFrame - centerX_nextFrame(1), interpolated_layer);
    %             centerY_nextFrame_interp = interpolateData(centerY_nextFrame - centerY_nextFrame(1), interpolated_layer);
    % 
    %             centerX_firstFrame_interp = centerX_firstFrame_interp(1:100:end);
    %             centerY_firstFrame_interp = centerY_firstFrame_interp(1:100:end);
    % 
    %             centerX_nextFrame_interp = centerX_nextFrame_interp(1:100:end);
    %             centerY_nextFrame_interp = centerY_nextFrame_interp(1:100:end);
    % 
    %             centerX_interpolationResult = centerX_gd(1) + centerX_firstFrame_interp*firstFrame_center_weight + centerX_nextFrame_interp*nextFrame_center_weight;
    %             centerY_interpolationResult = centerY_gd(1) + centerY_firstFrame_interp*firstFrame_center_weight + centerY_nextFrame_interp*nextFrame_center_weight;
    %     
    
    
end

