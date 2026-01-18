function thisEddyHistory=allpaths(subGraphMap,firstEddyNodeIndex,lastEddyNodeIndex)
%ALLPATHS Summary of this function goes here
%   Detailed explanation goes here

nodes = subGraphMap.Nodes;

data = nodes{firstEddyNodeIndex,"Name"};
splitData = strsplit(data{1},'.');
firstTimestep = str2double(splitData{1});

data = nodes{lastEddyNodeIndex,"Name"};
splitData = strsplit(data{1},'.');
lastTimestep = str2double(splitData{1});

for i = 1:1:size(nodes,1)
    data = nodes{i,"Name"};
    sequenceData = strsplit(data{1}, ".");
    eddyIndex(i,1) = str2double(sequenceData{1});
end

% thisEddyHistory = firstTimestep:1:lastTimestep;


% [~,firstLoc] = ismember(firstTimestep,eddyIndex);
% [~,lastLoc] = ismember(lastTimestep,eddyIndex);

eddyIndex = sort(eddyIndex);
thisEddyHistory = eddyIndex<=lastTimestep & eddyIndex>=firstTimestep;


% eddyIndex = sortcolumn(eddyIndex,1);

end