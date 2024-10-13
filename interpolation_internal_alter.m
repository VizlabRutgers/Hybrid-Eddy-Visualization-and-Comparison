function [tempEvaluationResult] = interpolation_internal_alter(wholecenterX, wholecenterY, wholeEddyDepth, interpolation_startIndex, originData_stride,z_val)
%INTERPOLATION_INTERNAL 此处显示有关此函数的摘要
%   此处显示详细说明

    
    
    centerX_firstFrame = wholecenterX{1, interpolation_startIndex};
    centerY_firstFrame = wholecenterY{1, interpolation_startIndex};
    radiusInDepth_firstFrame = wholeEddyDepth{1, interpolation_startIndex};
    
    centerX_nextFrame = wholecenterX{1, interpolation_startIndex+originData_stride};
    centerY_nextFrame = wholecenterY{1, interpolation_startIndex+originData_stride};
    radiusInDepth_nextFrame = wholeEddyDepth{1, interpolation_startIndex+originData_stride};
    
    tempEvaluationResult = [];

    % Decide the stride for interpolation
    % i.e. stride = 2 means one intermediate frame
    

    for interpolateIndex = 1:1:originData_stride-1

        centerX_gd = wholecenterX{1, interpolation_startIndex+interpolateIndex};
        centerY_gd = wholecenterY{1, interpolation_startIndex+interpolateIndex};
        radiusInDepth_gd = wholeEddyDepth{1, interpolation_startIndex+interpolateIndex};

        % Decide how many layers in depth
        interpolated_layer = length(radiusInDepth_firstFrame)+round((length(radiusInDepth_nextFrame)-length(radiusInDepth_firstFrame))/originData_stride*interpolateIndex);
        
        %interpolation of radiius
        
        if(length(radiusInDepth_firstFrame) == 1 || length(radiusInDepth_nextFrame) == 1)
            currEvaluationResult = {0,0,0,0,interpolation_startIndex, originData_stride,interpolation_startIndex+interpolateIndex,interpolated_layer,radiusInDepth_firstFrame(1),radiusInDepth_nextFrame(1),radiusInDepth_gd,0,0,0,0,0,0,0};
            tempEvaluationResult = [tempEvaluationResult;currEvaluationResult];
        else

            firstFrame_distance = sqrt((centerX_firstFrame(1)-centerX_gd(1))^2+(centerY_firstFrame(1)-centerY_gd(1))^2);
            nextFrame_distance = sqrt((centerX_nextFrame(1)-centerX_gd(1))^2+(centerY_nextFrame(1)-centerY_gd(1))^2);
    
            firstFrame_center_weight = nextFrame_distance / (firstFrame_distance + nextFrame_distance);
            nextFrame_center_weight = firstFrame_distance / (firstFrame_distance + nextFrame_distance);
    
            firstFrame_radius_weight = abs(radiusInDepth_gd(1) - radiusInDepth_nextFrame(1)) / (abs(radiusInDepth_gd(1) - radiusInDepth_firstFrame(1)) + abs(radiusInDepth_gd(1) - radiusInDepth_nextFrame(1)));
            nextFrame_radius_weight = abs(radiusInDepth_gd(1) - radiusInDepth_firstFrame(1)) / (abs(radiusInDepth_gd(1) - radiusInDepth_firstFrame(1)) + abs(radiusInDepth_gd(1) - radiusInDepth_nextFrame(1)));
    
            larger_layer = max(length(centerX_firstFrame), length(centerX_nextFrame));
    
            radius_firstFrame_interp = interp1(1:1:length(radiusInDepth_firstFrame), radiusInDepth_firstFrame,1:((length(radiusInDepth_firstFrame)-1)/(larger_layer-1)):length(radiusInDepth_firstFrame));
            radius_nextFrame_interp = interp1(1:1:length(radiusInDepth_nextFrame), radiusInDepth_nextFrame, 1:(length(radiusInDepth_nextFrame)-1)/(larger_layer-1):length(radiusInDepth_nextFrame));
            x_firstFrame_interp = interp1(1:1:length(radiusInDepth_firstFrame), centerX_firstFrame,1:((length(radiusInDepth_firstFrame)-1)/(larger_layer-1)):length(radiusInDepth_firstFrame));
            x_nextFrame_interp = interp1(1:1:length(radiusInDepth_nextFrame), centerX_nextFrame, 1:(length(radiusInDepth_nextFrame)-1)/(larger_layer-1):length(radiusInDepth_nextFrame));
            y_firstFrame_interp = interp1(1:1:length(radiusInDepth_firstFrame), centerY_firstFrame,1:((length(radiusInDepth_firstFrame)-1)/(larger_layer-1)):length(radiusInDepth_firstFrame));
            y_nextFrame_interp = interp1(1:1:length(radiusInDepth_nextFrame), centerY_nextFrame, 1:(length(radiusInDepth_nextFrame)-1)/(larger_layer-1):length(radiusInDepth_nextFrame));

    
    
    
            radius_depth_firstFrame = [1:((length(radiusInDepth_firstFrame)-1)/(larger_layer-1)):length(radiusInDepth_firstFrame);radius_firstFrame_interp];
            radius_depth_nextFrame = [1:(length(radiusInDepth_nextFrame)-1)/(larger_layer-1):length(radiusInDepth_nextFrame);radius_nextFrame_interp];
            x_depth_firstFrame = [1:((length(radiusInDepth_firstFrame)-1)/(larger_layer-1)):length(radiusInDepth_firstFrame);x_firstFrame_interp];
            x_depth_nextFrame = [1:((length(radiusInDepth_nextFrame)-1)/(larger_layer-1)):length(radiusInDepth_nextFrame);x_nextFrame_interp];
            y_depth_firstFrame = [1:(length(radiusInDepth_firstFrame)-1)/(larger_layer-1):length(radiusInDepth_firstFrame);y_firstFrame_interp];
            y_depth_nextFrame = [1:(length(radiusInDepth_nextFrame)-1)/(larger_layer-1):length(radiusInDepth_nextFrame);y_nextFrame_interp];




            radius_interpolationResult = radius_depth_firstFrame * (originData_stride-interpolateIndex) / originData_stride + radius_depth_nextFrame * interpolateIndex / originData_stride;
            x_interpolationResult = x_depth_firstFrame * (originData_stride-interpolateIndex) / originData_stride + x_depth_nextFrame * interpolateIndex / originData_stride;
            y_interpolationResult = y_depth_firstFrame * (originData_stride-interpolateIndex) / originData_stride + y_depth_nextFrame * interpolateIndex / originData_stride;



            radius_interpolationResult(2,:) = radius_interpolationResult(2,:) + (radiusInDepth_gd(1) - radius_interpolationResult(2,1));
            x_interpolationResult(2,:) = x_interpolationResult(2,:) + (centerX_gd(1) - x_interpolationResult(2,1));
            y_interpolationResult(2,:) = y_interpolationResult(2,:) + (centerY_gd(1) - y_interpolationResult(2,1));
    
            final_radius_interp = [];
            final_radius_interp(1,:) = 1:1:floor(max(radius_interpolationResult(1,:)));
            final_radius_interp(2,:) = interp1(radius_interpolationResult(1,:), radius_interpolationResult(2,:), 1:1:floor(max(radius_interpolationResult(1,:))));

            final_x_interp = [];
            final_y_interp = [];

            final_x_interp(1,:) = 1:1:floor(max(x_interpolationResult(1,:)));
            final_y_interp(1,:) = 1:1:floor(max(y_interpolationResult(1,:)));
            final_x_interp(2,:) = interp1(x_interpolationResult(1,:), x_interpolationResult(2,:), 1:1:floor(max(x_interpolationResult(1,:))));
            final_y_interp(2,:) = interp1(y_interpolationResult(1,:), y_interpolationResult(2,:), 1:1:floor(max(y_interpolationResult(1,:))));


    


       
    
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
    
    
%             figure,
%             p1 = plot(1:1:length(radiusInDepth_firstFrame),radiusInDepth_firstFrame,'r-o');
%             hold on;
%             p2 = plot(1:1:length(radiusInDepth_nextFrame),radiusInDepth_nextFrame,'b-o');
%             hold on;
%             p3 = plot(1:1:length(radiusInDepth_gd),radiusInDepth_gd,'k-o');
%             hold on;
%             p7 = plot( radius_interpolationResult(1,:), radius_interpolationResult(2,:),'k', "LineWidth",4);
%             hold on
%             p8 = plot( final_interp(1,:), final_interp(2,:),'g', "LineWidth",4);


%             hold off;
%             title("Ground truth for the eddy radius in three frames");
%             legend([p1, p2, p3, p7, p8], "frame"+int2str(interpolation_startIndex), "frame"+int2str(interpolation_startIndex+originData_stride), "frame"+int2str(interpolation_startIndex+interpolateIndex)+"(interpolation ground truth)", "fitting result", "fitting result on interger layer");
%             xlabel("depth");
%             ylabel("radius");
            dtwDistance = dtw(radiusInDepth_gd,final_radius_interp(2,:));
    
            [euclideanDistance,cosineDistance,corrDistance,ssimDistance] = measureSimilarity(final_radius_interp(2,:), radiusInDepth_gd');
            currEvaluationResult = {euclideanDistance,cosineDistance,corrDistance,ssimDistance,interpolation_startIndex, originData_stride,interpolation_startIndex+interpolateIndex,length(radiusInDepth_gd),radiusInDepth_firstFrame(1),radiusInDepth_nextFrame(1),radiusInDepth_gd,final_radius_interp(2,:)',dtwDistance, centerX_gd, centerY_gd, final_x_interp(2,:)',final_y_interp(2,:)'};

%             figure,
%             dtw(radiusInDepth_gd',final_interp(2,:));

    
            if(radiusInDepth_gd(1) < min(currEvaluationResult{9},currEvaluationResult{10}))
                currEvaluationResult{18} = abs(radiusInDepth_gd(1) - min(currEvaluationResult{9},currEvaluationResult{10}));
            elseif(radiusInDepth_gd(1) > max(currEvaluationResult{9},currEvaluationResult{10}))
                currEvaluationResult{18} = abs(radiusInDepth_gd(1) - max(currEvaluationResult{9},currEvaluationResult{10}));
            else
                currEvaluationResult{18} = 0;
            end


            tempEvaluationResult = [tempEvaluationResult;currEvaluationResult];

    %             figure,
    %             p1 = plot(1:1:length(radiusInDepth_firstFrame),radiusInDepth_firstFrame,'r-o');
    %             hold on;
    %             p2 = plot(1:1:length(radiusInDepth_nextFrame),radiusInDepth_nextFrame,'b-o');
    %             hold on;
    %             p3 = plot(1:0.01:interpolated_layer,radius_firstFrame_interp,'r--');
    %             hold on;
    %             p4 = plot(1:0.01:interpolated_layer,radius_nextFrame_interp,'b--');
    %             hold on;
    %             p5 = plot(1:0.01:interpolated_layer,radius_firstFrame_interp_adjusted,'r:.');
    %             hold on;
    %             p6 = plot(1:0.01:interpolated_layer,radius_nextFrame_interp_adjusted,'b:.');
    %             hold on;
    %             p7 = plot(1:0.01:interpolated_layer,radius_interpolationResult,'k:.');
    %             hold on;
    %     
    %             p8 = plot(1:1:length(radiusInDepth_gd),radiusInDepth_gd,'k-o');
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
    %         p2 = plot(1:1:length(radiusInDepth_nextFrame),radiusInDepth_nextFrame,'b-o');
    %         hold on;
    % %         p3 = plot(1:0.01:interpolated_layer,radius_firstFrame_interp,'r--');
    % %         hold on;
    % %         p4 = plot(1:0.01:interpolated_layer,radius_nextFrame_interp,'b--');
    % %         hold on;
    %         p5 = plot(1:0.01:interpolated_layer,radius_firstFrame_interp_adjusted,'r:.');
    %         hold on;
    %         p6 = plot(1:0.01:interpolated_layer,radius_nextFrame_interp_adjusted,'b:.');
    %         hold on;
    %         p7 = plot(1:0.01:interpolated_layer,radius_interpolationResult,'k:.');
    %         hold on;
    % 
    %         p8 = plot(1:1:length(radiusInDepth_gd),radiusInDepth_gd,'k-o');
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
    %     %         p2 = plot(1:1:length(radiusInDepth_nextFrame),radiusInDepth_nextFrame,'b-o');
    %     %         hold on;
    %     %         p3 = plot(1:0.01:interpolated_layer,radius_firstFrame_interp,'r--');
    %     %         hold on;
    %     %         p4 = plot(1:0.01:interpolated_layer,radius_nextFrame_interp,'b--');
    %     %         hold on;
    %     %         p5 = plot(1:0.01:interpolated_layer,radius_firstFrame_interp_adjusted,'r:.');
    %     %         hold on;
    %     %         p6 = plot(1:0.01:interpolated_layer,radius_nextFrame_interp_adjusted,'b:.');
    %     %         hold on;
    %             p7 = plot(1:0.01:interpolated_layer,radius_interpolationResult,'k:.');
    %             hold on;
    %     
    %             p8 = plot(1:1:length(radiusInDepth_gd),radiusInDepth_gd,'k-o');
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

