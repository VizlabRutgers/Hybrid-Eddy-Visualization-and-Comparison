function [vq1] = interpolateData(radiusInDepth_firstFrame,interpolated_layer)
    
    x = 1:1:length(radiusInDepth_firstFrame);
    x_expanded = 1:(interpolated_layer-1) / (length(radiusInDepth_firstFrame)-1):interpolated_layer;
    v = radiusInDepth_firstFrame;
    xq = 1:0.01:interpolated_layer;
    
    vq1 = interp1(x_expanded,v*interpolated_layer/length(radiusInDepth_firstFrame),xq);

%     plot(x,v,':o',xq,vq1,':.');
%     xlim([0 2*pi]);
%     title('(Default) Linear Interpolation');
end

