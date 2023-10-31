function [eddyHistory] = eddyPathProcess(clockEddy,conterclockEddy,eddyPath,eddyIndex,G,SizeLimits)
%clockEddy: clockwise eddy data
% conterclockEddy:counterclockwise eddy data
% eddyPath: eddy's paths from Graph
% eddyIndex: The eddy index defined by user
% G: eddy graph
% SizeLimits: minimum of eddy size. eddy below this size will be filtered out 
%   Detailed explanation goes here
    clcIndex=1;
    counterIndex=1;
    NotRepeatFlag=0;
    eddyHistory = {};



    

    for i = 1:1:length(eddyIndex)
        thisEddyHistoryIndex=eddyIndex(i);

        subpath=eddyPath{thisEddyHistoryIndex};
        subGraphMap=subgraph(G,subpath);
        eddyIndividualPaths=getpaths(subGraphMap);
        subEddyNodes=table2cell(subGraphMap.Nodes);


%         figure,
%         plot(subGraphMap)
%         title(['structure of eddypath ' num2str(thisEddyHistoryIndex)])
%         ylabel('time progresses downwards')
%         xlabel('node linages')


        for pathIndex=1:1:length(eddyIndividualPaths)
            individualPath=eddyIndividualPaths{pathIndex};
            individualPath=arrayfun(@(x) subEddyNodes(x), individualPath)';
            if(size(individualPath,1)<=SizeLimits)
                continue;
            end
            thisEddy=cellfun(@(x) strsplit(x, '.'), individualPath, 'UniformOutput', false);
            thisEddy=str2double(vertcat(thisEddy{:})).';

            eddyHistoryTemp=[];


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
                    error('Error: Not clockwise or conterclockwise');
                end

            end

            

        end
        eddyHistory{end+1}=eddyHistoryTemp;
        
    end
    eddyHistory=eddyHistory';
end

