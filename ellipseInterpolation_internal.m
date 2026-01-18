function [InterpResult] = ellipseInterpolation_internal(wholecenterX, wholecenterY, wholeLengthX, wholeLengthY,wholeTheta, interpolation_startIndex, originData_stride,z_val,eddyPathNumber)
%INTERPOLATION_INTERNAL 此处显示有关此函数的摘要
%   此处显示详细说明

    
    
    centerX_previousFrame = wholecenterX{1, interpolation_startIndex};
    centerY_previousFrame = wholecenterY{1, interpolation_startIndex};
    lengthX_previousFrame = wholeLengthX{1, interpolation_startIndex};
    lengthY_previousFrame = wholeLengthY{1, interpolation_startIndex};
    theta_previousFrame = wholeTheta{1, interpolation_startIndex};
    
    centerX_nextFrame = wholecenterX{1, interpolation_startIndex+originData_stride};
    centerY_nextFrame = wholecenterY{1, interpolation_startIndex+originData_stride};
    lengthX_nextFrame = wholeLengthX{1, interpolation_startIndex+originData_stride};
    lengthY_nextFrame = wholeLengthY{1, interpolation_startIndex+originData_stride};
    theta_nextFrame = wholeTheta{1, interpolation_startIndex+originData_stride};
    
    InterpResult = {};

    % Decide the stride for interpolation
    % i.e. stride = 2 means one intermediate frame
    


    for interpolateIndex = 1:1:originData_stride-1

        centerX_gd = wholecenterX{1, interpolation_startIndex+interpolateIndex};
        centerY_gd = wholecenterY{1, interpolation_startIndex+interpolateIndex};
        lengthX_gd = wholeLengthX{1, interpolation_startIndex+interpolateIndex};
        lengthY_gd = wholeLengthY{1, interpolation_startIndex+interpolateIndex};
        theta_gd = wholeTheta{1, interpolation_startIndex+interpolateIndex};

        % Decide how many layers in depth
        interpolated_layer = length(lengthX_previousFrame)+round((length(lengthX_nextFrame)-length(lengthX_previousFrame))/originData_stride*interpolateIndex);
        
        %% interpolation of radiius
        
        % exclude the case that only has one layer (on the surface)
        if(length(lengthX_previousFrame) == 1 || length(lengthX_nextFrame) == 1)
            currEvaluationResult = {0,0,0,0,interpolation_startIndex, originData_stride,interpolation_startIndex+interpolateIndex,interpolated_layer,lengthX_previousFrame(1),lengthX_nextFrame(1),lengthX_gd,0,0,0,0,0,0,0};
            evaluationResult = [evaluationResult;currEvaluationResult];
        else

            % calculate the distance between middleFrame-previousFrame &
            % middleFrame-nextFrame
%             firstFrame_distance = sqrt((centerX_previousFrame(1)-centerX_gd(1))^2+(centerY_previousFrame(1)-centerY_gd(1))^2);
%             nextFrame_distance = sqrt((centerX_nextFrame(1)-centerX_gd(1))^2+(centerY_nextFrame(1)-centerY_gd(1))^2);
%     
%             firstFrame_center_weight = nextFrame_distance / (firstFrame_distance + nextFrame_distance);
%             nextFrame_center_weight = firstFrame_distance / (firstFrame_distance + nextFrame_distance);
%     
%             firstFrame_lengthX_weight = abs(lengthX_gd(1) - lengthX_nextFrame(1)) / (abs(lengthX_gd(1) - lengthX_previousFrame(1)) + abs(lengthX_gd(1) - lengthX_nextFrame(1)));
%             nextFrame_lengthX_weight = abs(lengthX_gd(1) - lengthX_previousFrame(1)) / (abs(lengthX_gd(1) - lengthX_previousFrame(1)) + abs(lengthX_gd(1) - lengthX_nextFrame(1)));
%     
            % The deepest layer between two frames
   
            deeperLayer = max(length(centerX_previousFrame), length(centerX_nextFrame));
            depthPrev = 1:((length(centerX_previousFrame)-1)/(deeperLayer-1)):length(centerX_previousFrame);
            depthNext = 1:(length(centerX_nextFrame)-1)/(deeperLayer-1):length(centerX_nextFrame);
            w1 = (originData_stride - interpolateIndex) / originData_stride;
            w2 = interpolateIndex / originData_stride;
            interp_depth = depthPrev * w1 + depthNext * w2;

            final_lengthX_interp = interpFunc(lengthX_previousFrame, lengthX_nextFrame, lengthX_gd,originData_stride,interpolateIndex,interp_depth);
            final_lengthY_interp = interpFunc(lengthY_previousFrame, lengthY_nextFrame, lengthY_gd,originData_stride,interpolateIndex,interp_depth);
            final_theta_interp = interpAngleFunc(theta_previousFrame, theta_nextFrame, theta_gd,originData_stride,interpolateIndex,interp_depth);
            final_centerX_interp = interpFunc(centerX_previousFrame, centerX_nextFrame, centerX_gd,originData_stride,interpolateIndex,interp_depth);
            final_centerY_interp = interpFunc(centerY_previousFrame, centerY_nextFrame, centerY_gd,originData_stride,interpolateIndex,interp_depth);
    


       
    

%             figure,
%             p1 = plot(1:1:length(radiusInDepth_firstFrame),radiusInDepth_firstFrame,'r-o');
%             hold on;
%             p2 = plot(1:1:length(lengthX_nextFrame),lengthX_nextFrame,'b-o');
%             hold on;
%             p3 = plot(1:1:length(lengthX_gd),lengthX_gd,'k-o');
%             hold on;
%             p7 = plot( lengthX_interpolationResult(1,:), lengthX_interpolationResult(2,:),'k', "LineWidth",4);
%             hold on
%             p8 = plot( final_interp(1,:), final_interp(2,:),'g', "LineWidth",4);


%             hold off;
%             title("Ground truth for the eddy radius in three frames");
%             legend([p1, p2, p3, p7, p8], "frame"+int2str(interpolation_startIndex), "frame"+int2str(interpolation_startIndex+originData_stride), "frame"+int2str(interpolation_startIndex+interpolateIndex)+"(interpolation ground truth)", "fitting result", "fitting result on interger layer");
%             xlabel("depth");
%             ylabel("radius");

            currEvaluationResult = ellipseEvaluation(centerX_gd, centerY_gd,lengthX_gd, lengthY_gd, theta_gd, ...
                final_centerX_interp', final_centerY_interp',final_lengthX_interp', final_lengthY_interp', final_theta_interp');
            
            gd_data = struct();
            gd_data.centerX = centerX_gd;
            gd_data.centerY = centerY_gd;
            gd_data.lengthX = lengthX_gd;
            gd_data.lengthY = lengthY_gd;
            gd_data.theta = theta_gd;

            interp_data = struct();
            interp_data.centerX = final_centerX_interp';
            interp_data.centerY = final_centerY_interp';
            interp_data.lengthX = final_lengthX_interp';
            interp_data.lengthY = final_lengthY_interp';
            interp_data.theta = final_theta_interp';

            currInterpResult = struct();
            currInterpResult.startIndex = interpolation_startIndex;
            currInterpResult.interpolationStride = originData_stride;
            currInterpResult.currentInterpolationIndex = interpolation_startIndex+interpolateIndex;
            currInterpResult.endIndex = interpolation_startIndex+originData_stride;
            currInterpResult.gd_data = gd_data; 
            currInterpResult.interp_data = interp_data; 
            currInterpResult.evaluationMetric = currEvaluationResult;
            currInterpResult.eddyPathNumber = eddyPathNumber;
            
            
            InterpResult = [InterpResult;currInterpResult];

    %             figure,
    %             p1 = plot(1:1:length(radiusInDepth_firstFrame),radiusInDepth_firstFrame,'r-o');
    %             hold on;
    %             p2 = plot(1:1:length(lengthX_nextFrame),lengthX_nextFrame,'b-o');
    %             hold on;
    %             p3 = plot(1:0.01:interpolated_layer,lengthX_firstFrame_interp,'r--');
    %             hold on;
    %             p4 = plot(1:0.01:interpolated_layer,lengthX_nextFrame_interp,'b--');
    %             hold on;
    %             p5 = plot(1:0.01:interpolated_layer,lengthX_firstFrame_interp_adjusted,'r:.');
    %             hold on;
    %             p6 = plot(1:0.01:interpolated_layer,lengthX_nextFrame_interp_adjusted,'b:.');
    %             hold on;
    %             p7 = plot(1:0.01:interpolated_layer,lengthX_interpolationResult,'k:.');
    %             hold on;
    %     
    %             p8 = plot(1:1:length(lengthX_gd),lengthX_gd,'k-o');
    %     
    %             hold off;
    %             title("Interpolation for the eddy radius in three frames");
    %             legend([p1, p2, p8, p3, p4, p5, p6,p7], ...
    %                 "first frame ground truth", "third frame ground truth", "second frame ground truth",...
    %                 "first frame interpolation", "third frame interpolation", ...
    %                 "first frame interpolation after adjusting", "third frame interpolation after adjusting",...
    %                 "fitting result");
    %     
    
    
    %         figure,
    %         p1 = plot(1:1:length(radiusInDepth_firstFrame),radiusInDepth_firstFrame,'r-o');
    %         hold on;
    %         p2 = plot(1:1:length(lengthX_nextFrame),lengthX_nextFrame,'b-o');
    %         hold on;
    % %         p3 = plot(1:0.01:interpolated_layer,lengthX_firstFrame_interp,'r--');
    % %         hold on;
    % %         p4 = plot(1:0.01:interpolated_layer,lengthX_nextFrame_interp,'b--');
    % %         hold on;
    %         p5 = plot(1:0.01:interpolated_layer,lengthX_firstFrame_interp_adjusted,'r:.');
    %         hold on;
    %         p6 = plot(1:0.01:interpolated_layer,lengthX_nextFrame_interp_adjusted,'b:.');
    %         hold on;
    %         p7 = plot(1:0.01:interpolated_layer,lengthX_interpolationResult,'k:.');
    %         hold on;
    % 
    %         p8 = plot(1:1:length(lengthX_gd),lengthX_gd,'k-o');
    % 
    %         hold off;
    %         title("Interpolation for the eddy radius in three frames");
    %         legend([p1, p2, p8, p5, p6,p7], ...
    %             "first frame ground truth", "third frame ground truth", "second frame ground truth",...
    %             "first frame interpolation after adjusting", "third frame interpolation after adjusting",...
    %             "fitting result");
    
    
    
    
    
    %             figure,
    %     %         p1 = plot(1:1:length(radiusInDepth_firstFrame),radiusInDepth_firstFrame,'r-o');
    %     %         hold on;
    %     %         p2 = plot(1:1:length(lengthX_nextFrame),lengthX_nextFrame,'b-o');
    %     %         hold on;
    %     %         p3 = plot(1:0.01:interpolated_layer,lengthX_firstFrame_interp,'r--');
    %     %         hold on;
    %     %         p4 = plot(1:0.01:interpolated_layer,lengthX_nextFrame_interp,'b--');
    %     %         hold on;
    %     %         p5 = plot(1:0.01:interpolated_layer,lengthX_firstFrame_interp_adjusted,'r:.');
    %     %         hold on;
    %     %         p6 = plot(1:0.01:interpolated_layer,lengthX_nextFrame_interp_adjusted,'b:.');
    %     %         hold on;
    %             p7 = plot(1:0.01:interpolated_layer,lengthX_interpolationResult,'k:.');
    %             hold on;
    %     
    %             p8 = plot(1:1:length(lengthX_gd),lengthX_gd,'k-o');
    %     
    %             hold off;
    %             title("Interpolation for the eddy radius in three frames");
    %             legend([p8, p7], ...
    %                 "second frame ground truth", "fitting result");
    % 
    %             figure,
    %             plot3(centerX_firstFrame, centerY_firstFrame, z_val(1:length(centerX_firstFrame)),'r-o');
    %             hold on;
    %             plot3(centerX_nextFrame, centerY_nextFrame, z_val(1:length(centerX_nextFrame)),'b-o');
    %             hold on;
    %             plot3(centerX_gd, centerY_gd, z_val(1:length(centerX_gd)),'k-o');
    %             hold on;
    %             plot3(centerX_interpolationResult, centerY_interpolationResult, z_val(1:length(centerX_interpolationResult)),'k-', "LineWidth",5);
    %             daspect([1 1 450]);
    %             xlabel('Longitude');
    %             ylabel('Latitude');
    %             zlabel("Depth (meter)")
    % 
    %             set(gca,"ZDir", "reverse");
    

        end
    end
end

