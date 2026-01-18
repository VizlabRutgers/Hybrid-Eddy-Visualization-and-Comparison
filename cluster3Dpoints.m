function [finalTracks, eddy_ellipse_objects] = cluster3Dpoints(points, distanceThreshold,track_depthThres,x_val,y_val,z_val,dataResolution)
% points: N x 3 matrix, each row is (x, y, z)
% distanceThreshold: max distance to associate points across layers
% tracks: a cell array, each cell is a list of points for one target


% 提取所有 z 值并排序
z_values = 1:1:length(z_val);
num_layers = length(points);

% 初始化轨迹，第一层每个点一条轨迹
finalTracks = {};
activeTracks = {};
surfacePoints = points{z_values(1)};

surfacePoints = cellfun(@(x) num2cell(x), surfacePoints, 'UniformOutput', false);
surfacePoints  = vertcat(surfacePoints{:});   

previousTrackPtMap = [];

for i = 1:size(surfacePoints,2)
    activeTracks{end+1} = surfacePoints(:,i);
    previousTrackPtMap(i,:) = cell2mat(surfacePoints(1:2,i))';
end

% 从第二层开始逐层匹配
for i = 2:num_layers

    z = z_values(i);
    current_pts = points{i};
    if numel(current_pts{1})==0 
        finalTracks = [finalTracks,activeTracks];
        break;
    end
    
    if isempty(activeTracks)
        break;
    end

    current_pts = cellfun(@(x) num2cell(x), current_pts, 'UniformOutput', false);
    current_pts  = vertcat(current_pts{:}); 
    
    if(i==7)
        i=7;
    end
    num_pts = size(current_pts, 2);
    pt_assigned = false(num_pts, 1);
    track_assigned = false(length(activeTracks),1);
    
    currentTrackPtMap = cell2mat(current_pts(1:2,:))';
    M = size(currentTrackPtMap,1);
    N = size(previousTrackPtMap,1);
    
    previousTrackPtMap = cellfun(@(x) x(1:2,end), activeTracks, 'UniformOutput',false);
    previousTrackPtMap = cell2mat(cellfun(@(x) [x{1},x{2}], previousTrackPtMap,'UniformOutput',false)');

    distanceMap = pdist2(currentTrackPtMap, previousTrackPtMap, 'euclidean');
    
    
    
    % 
    % L = max(M, N);
    % cost = ones(L, L) * (distanceThreshold*10);  % 大值填充
    % cost(1:M, 1:N) = distanceMap;
    
    % 4. 匈牙利算法求解最小分配
    % MATLAB 自带：matchpairs （需要 Statistics and Machine Learning Toolbox）
    [matches,unMatchCurrent,unMatchPrevious] = matchpairs(distanceMap, distanceThreshold);
    


    for matchIndex = 1:size(matches, 1)
        previousTrack_temp = activeTracks{matches(matchIndex,2)};
        currentPt_temp = current_pts(:,matches(matchIndex,1));


        activeTracks{matches(matchIndex,2)} = [previousTrack_temp, currentPt_temp];
    end
        
    if(~isempty(unMatchPrevious))
        for unMatchIndex = 1:1:length(unMatchPrevious)
            unMatchTrack = activeTracks{unMatchPrevious(unMatchIndex)};
            if(size(unMatchTrack, 2) >= track_depthThres)
                finalTracks(end+1) = {unMatchTrack};
            end
        end
        activeTracks(unMatchPrevious) = [];
    end

    
%     if(~isempty(unMatchCurrent))
%      
%         newTrack = mat2cell(current_pts(:, unMatchCurrent)', ...
%             repmat([1],[length(unMatchCurrent),1]));
%         newTrack = cellfun(@transpose,newTrack, "UniformOutput",false);
%         activeTracks = [activeTracks, newTrack'];
%     end

    % % Not assigned track is regarded as "end" at current layer.
    % for t = 1:length(activeTracks)
    %     if (~track_assigned(t))
    %         finalTracks{end+1} = activeTracks{t};
    %         activeTracks{t} = [];
    %     end
    % end
    % activeTracks = activeTracks(~cellfun(@isempty, activeTracks));
    % 
    % % 把没有被匹配的点新建成轨迹（可能是新目标或噪声）
    % for j = 1:num_pts
    %     if ~pt_assigned(j)
    %         activeTracks{end+1} = current_pts(:,j);
    %     end
    % end




    % % 为每条轨迹找最近的点
    % for t = 1:length(activeTracks)
    %     currentTrack = activeTracks{t};
    %     last_pt = currentTrack(:,end);
    %     min_dist = Inf;
    %     best_idx = -1;
    % 
    %     for j = 1:num_pts
    %         if pt_assigned(j)
    %             continue;
    %         end
    %         d = norm(current_pts(1:2,j) - last_pt(1:2));
    %         if d < distanceThreshold && d < min_dist
    %             min_dist = d;
    %             best_idx = j;
    %         end
    %     end
    % 
    %     if best_idx ~= -1
    %         activeTracks{t} = [activeTracks{t}, current_pts(:,best_idx)];
    %         pt_assigned(best_idx) = true;
    %         track_assigned(t) = true;
    %     end
    % end




    % % Not assigned track is regarded as "end" at current layer.
    % for t = 1:length(activeTracks)
    %     if (~track_assigned(t))
    %         finalTracks{end+1} = activeTracks{t};
    %         activeTracks{t} = [];
    %     end
    % end
    % activeTracks = activeTracks(~cellfun(@isempty, activeTracks));
    % 
    % % 把没有被匹配的点新建成轨迹（可能是新目标或噪声）
    % for j = 1:num_pts
    %     if ~pt_assigned(j)
    %         activeTracks{end+1} = current_pts(:,j);
    %     end
    % end
    % 
    % 
    % if i == num_layers
    %     finalTracks = [finalTactiveTracks,activeTracks];
    % end


end




%% Gather the pts for each individual eddy
% close all;

eddy_ellipse_objects = {};

% fh1 = figure();
% fh1.WindowState = 'maximized';
% ax1=axes(fh1);
% set(gca,'FontSize',16);

% coastFaces = load('coastFaceData_RS.mat','coastFaces').coastFaces;
% coastVerts = load('coastVertData_RS.mat','coastVerts').coastVerts;
% seabed = patch('Faces',coastFaces, 'Vertices',coastVerts, ...  
% 'FaceColor',[.8, .8, .8],...
% 'EdgeColor','none',...
% 'FaceAlpha',0.5, ...
% 'Parent', ax1);
% hold on



% for trackIndex=26
for trackIndex=1:1:length(finalTracks)
    eddy_ellipse_pts = [];
    currentTrack = finalTracks{trackIndex};

    for depthIndex = 1:1:size(currentTrack,2)

        

        
        eddy_center_x = currentTrack{1,depthIndex};
        eddy_center_y = currentTrack{2,depthIndex};
        
        [x_grid,y_grid] = meshgrid(x_val,y_val);
%         if(eddy_length_x < dataResolution && eddy_length_y <dataResolution)
%             break;
%         end   
        
        if(size(currentTrack,1)==5)
            ellipse_x_shift = currentTrack{4,depthIndex}{1};
            ellipse_y_shift = currentTrack{5,depthIndex}{1};
            eddy_elipse_x = eddy_center_x + ellipse_x_shift;
            eddy_elipse_y = eddy_center_y + ellipse_y_shift;
%             plot3(ax1,eddy_elipse_x,eddy_elipse_y,z_val(depthIndex*ones(360,1)),'r.');
            [in,~] = inpolygon(x_grid,y_grid,eddy_elipse_x,eddy_elipse_y);
            filledElipsoidPoints_X = x_grid(in);
            filledElipsoidPoints_Y = y_grid(in);
            
        
        
        elseif(size(currentTrack,1)==7)
            eddy_elipse_x = currentTrack{6,depthIndex}{1};
            eddy_elipse_y = currentTrack{7,depthIndex}{1};
%             plot3(ax1,eddy_elipse_x,eddy_elipse_y,z_val(depthIndex*ones(360,1)),'r.');
            [in,~] = inpolygon(x_grid,y_grid,eddy_elipse_x,eddy_elipse_y);
            filledElipsoidPoints_X = x_grid(in);
            filledElipsoidPoints_Y = y_grid(in);     
            
        else
            eddyInfoOnLayer = currentTrack(:, depthIndex);
            eddy_ellipse_theta =eddyInfoOnLayer(7);
            eddy_length_x = eddyInfoOnLayer(5);
            eddy_length_y = eddyInfoOnLayer(6);



            %     d(:,1);
            %     d(:,2);
            it = 1:360;
            ellipse_x_shift = eddy_length_x*cosd(it);
            ellipse_y_shift = eddy_length_y*sind(it);
            % Create rotation matrix
            R = [cosd(eddy_ellipse_theta) -sind(eddy_ellipse_theta); sind(eddy_ellipse_theta) cosd(eddy_ellipse_theta)];
            % Rotate points
            eddypts = R*[ellipse_x_shift ; ellipse_y_shift];
            ellipse_x_shift = eddypts(1,:);
            ellipse_y_shift = eddypts(2,:);
            eddy_elipse_x = eddy_center_x + ellipse_x_shift;
            eddy_elipse_y = eddy_center_y + ellipse_y_shift;

%             plot3(ax1,eddy_elipse_x,eddy_elipse_y,z_val(depthIndex*ones(360,1)),'r.')

            [filledElipsoidPoints_X,filledElipsoidPoints_Y] = fillElipsoidInsidePoints(eddy_center_x, eddy_center_y,...,
                eddy_elipse_x, eddy_elipse_y, eddy_length_x, eddy_length_y, eddy_ellipse_theta,x_val,y_val);
        end
        dataSize = size(filledElipsoidPoints_Y,1);

        if(~isempty(filledElipsoidPoints_X))

%         scatter(filledElipsoidPoints_X, filledElipsoidPoints_Y,"b.");
            eddy_ellipse_dataOnLayer = [eddy_center_x*ones(dataSize,1), eddy_center_y*ones(dataSize,1), filledElipsoidPoints_X,filledElipsoidPoints_Y,depthIndex*ones(dataSize,1)];
            eddy_ellipse_pts = [eddy_ellipse_pts;eddy_ellipse_dataOnLayer];

        end
        
    end 
    
    
%     if(~isempty(eddy_ellipse_pts))
%         scatter3(ax1,eddy_ellipse_pts(:,3),eddy_ellipse_pts(:,4),z_val(eddy_ellipse_pts(:,5)), "b.");
%         view(3);
%         hold on   
%     end

    eddy_ellipse_objects = [eddy_ellipse_objects, eddy_ellipse_pts];
end
% view(3);
% set(ax1, "ZDir", "reverse");
% zlim([0,500]);

end
