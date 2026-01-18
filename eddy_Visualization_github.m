clc;
clear all;
% close all;
% 

% This program is used to process the eddy information after the feature
% tracking program. The eddyVis and eddyGlobalVis function could also get
% the visualization and statistics.

%% The path of Feature tracking result

% The path of source dataset file (.nc file)
srcData="../../redsea/0001/COMBINED_2011013100.nc";

% Output from the Hybrid Eddy Detection
dataFilePath = "../../FT_result/1/"; 
trackTablePath = dataFilePath+dir(fullfile(dataFilePath,"*.trakTable")).name;

% The path of ow detection file (from Feature Tracking detection)
owDataFilePath="../../FT_OW/1/";

% The path of windingAngle detection file (from windingAngle detection)
windingAngleDataFilePath="./red_sea_windingAngle_fram_1.mat";


%% Program starts here
% ----------------------------------------------
%----------------------------------------------

allEddy=[];
clockEddy = [];
conterclockEddy = [];

%% Get eddy detection information from Feature tracking code
% The variable "allEddy" includes all the eddy statistic information

% clockwise eddy structures
Path = dataFilePath+"Seperated Structures/clockwise/";
file = dir(fullfile(Path,'*.uocd'));
filenames = {file.name};

for fileIndex = 1:1:length(filenames)
    ReadFileName = Path+filenames{fileIndex};
    data = load(ReadFileName);
    
    data(:,17)=(data(:,12)-3)./(data(:,13)-3);

%     radiusRange = unique(data(:,17));
%     meandata=[];
%     for i=1:1:length(radiusRange)
%         radiusData = data(data(:,17)==radiusRange(i),:);
%         % salt：10, temp：11
%         meandata(i)=mean(radiusData(:,11));
%     end
%     plot(radiusRange,meandata,"k","linewidth",2);
%     hold on

    clockEddy=[clockEddy;data(1,:)];
    

    allEddy=[allEddy;data(1,:)];
end

% counter-clockwise eddy structures
Path = dataFilePath+"Seperated Structures/counterclockwise/";
file = dir(fullfile(Path,'*.uocd'));
filenames = {file.name};

for fileIndex = 1:1:length(filenames)
    ReadFileName = Path+filenames{fileIndex};
    data = load(ReadFileName);

    data(:,17)=(data(:,12)-3)./(data(:,13)-3);

    radiusRange = unique(data(:,17));
%     meanSalt=[];
%     for i=1:1:length(radiusRange)
%         radiusData = data(data(:,17)==radiusRange(i),:);
%         meanSalt(i)=mean(radiusData(:,10));
%     end
%     plot(radiusRange,meanSalt,"k","linewidth",2);
%     hold on

    conterclockEddy=[conterclockEddy;data(1,:)];
    allEddy=[allEddy;data(1,:)];
end


%----------------------------------------------
%----------------------------------------------

%% Compute the eddy tracking information
[eddyPath,eddyNodes,G] = eddyActivity(trackTablePath);

%% Get individual eddy history information 
% Get the eddy history variable for an individual eddy. 
% for i = 1:1:size(eddyPath,1)
%     thisPath=eddyPath(i);
%     FindAEddy(i)=cellfun(@(y) contains("1.17",y),thisPath);
% end
% eddyExists=find(FindAEddy==1);

% eddyPathNum=118;
eddyPathNum=1:1:length(eddyPath);
% eddyPathNum=11;
eddyHistory = eddyPathProcess(clockEddy,conterclockEddy,eddyPath,eddyPathNum,G,1);



%% Visualize single path
% You could use this to visualize an individual eddy
% You need to change the variable name for the NC file inside the function

eddyVis(eddyPathNum,G,eddyPath,eddyHistory,dataFilePath,srcData);

%% Eddy global comparison
frameIndex = 1;
eddyMethodComparison(owDataFilePath,windingAngleDataFilePath,srcData, allEddy, frameIndex);
