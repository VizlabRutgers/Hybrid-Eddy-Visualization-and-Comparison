function [eddyHistory] = eddyPathProcess(clockEddy,conterclockEddy,eddyPath, ...
    G,livingFrameLimit)
%clockEddy: clockwise eddy data
% conterclockEddy:counterclockwise eddy data
% eddyPath: eddy's paths from Graph
% eddyIndex: The eddy index defined by user
% G: eddy graph
% SizeLimits: minimum of eddy size. eddy below this size will be filtered out 
%   Detailed explanation goes here

    eddyHistory = {};



    

    for i = 1:1:length(eddyPath)
        if(i==180)
            test = 0;
        end

        subpath=eddyPath{i};
        subGraphMap=subgraph(G,subpath);
        eddyIndividualPaths=getpaths(subGraphMap);
        subEddyNodes=table2cell(subGraphMap.Nodes);
        eddyHistoryTemp=[];

%         figure,
%         plot(subGraphMap)
%         title(['structure of eddypath ' num2str(thisEddyHistoryIndex)])
%         ylabel('time progresses downwards')
%         xlabel('node linages')


        for pathIndex=1:1:length(eddyIndividualPaths)
            individualPath=eddyIndividualPaths{pathIndex};
            individualPath=arrayfun(@(x) subEddyNodes(x), individualPath)';
            if(size(individualPath,1)<=livingFrameLimit)
                continue;
            end

            thisEddy=cellfun(@(x) strsplit(x, '.'), individualPath, 'UniformOutput', false);
            thisEddy=str2double(vertcat(thisEddy{:})).';




            [~, eddyTimeFrames] = size(thisEddy);
            for thisEddyIndex=1:1:eddyTimeFrames
                % Eddy rotation check
                clockwiseCheck=clockEddy((clockEddy(:,15)==thisEddy(1,thisEddyIndex))&(clockEddy(:,16)==thisEddy(2,thisEddyIndex)-1),:);
                conclockwiseCheck=conterclockEddy((conterclockEddy(:,15)==thisEddy(1,thisEddyIndex))&(conterclockEddy(:,16)==thisEddy(2,thisEddyIndex)-1),:);


                if(~isempty(clockwiseCheck)&&isempty(conclockwiseCheck))
                    eddyHistoryTemp = [eddyHistoryTemp;clockwiseCheck];
                elseif(~isempty(conclockwiseCheck)&&isempty(clockwiseCheck))
                    eddyHistoryTemp = [eddyHistoryTemp;conclockwiseCheck];
                else
                    warning("Something wrong in eddy history")
                end

            end

            

        end
        if(~isempty(eddyHistoryTemp))
            eddyHistory{end+1}=eddyHistoryTemp;
        end
        
    end
    eddyHistory=eddyHistory';
end

