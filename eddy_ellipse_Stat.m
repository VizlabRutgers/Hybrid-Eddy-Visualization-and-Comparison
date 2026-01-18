clc;
clear all;
% close all;
% 

% This program is used to process the eddy information after the feature
% tracking program. The eddyVis and eddyGlobalVis function could also get
% the visualization and statistics.

%% The path of Feature tracking result


% dataFilePath = "../../NA_5Day_result/";
% Path = dataFilePath+"Seperated Structures/clockwise/";
% trackTablePath = dataFilePath+dir(fullfile(dataFilePath,'*.trakTable')).name;
% % The path of source file (.nc file)
% srcData="../../NA_dataset/19980701.ocean_5day.nc";


% dataFilePath = "../Dataset/North_Pacific/result/";
% Path = dataFilePath+"Seperated Structures/clockwise/";
% trackTablePath = dataFilePath+dir(fullfile(dataFilePath,'*.trakTable')).name;
% % The path of source file (.nc file)
% srcData="../Dataset/North_Pacific/GLBy0_expt_93.nc4";

% dataFilePath = "../../NA_dataset/NA_new_1_25_cropped/19960101_cropped/"; % Output from the Hybrid Eddy Detection
% Path = dataFilePath+"Seperated Structures/clockwise/";
% trackTablePath = dataFilePath+dir(fullfile(dataFilePath,"*.trakTable")).name;
% % The path of source file (.nc file)
% srcData="../../NA_dataset/19960101.ocean_5day_cropped.nc";
% owDataFilePath="../../FT_OW/1/";


%% MexicoBay
% dataFilePath = "../../MexicoBay/windingAngle/FT_result_oldWindingAngle/"; % Output from the windingAngle Eddy Detection
% dataFilePath_new = "../../MexicoBay/windingAngle/FT_result_newWindingAngle/"; % Output from the windingAngle Eddy Detection
% dataFilePath_hybrid = "../../MexicoBay/hybrid/FT_result/";
% 
% Path = dataFilePath+"Seperated Structures/";
% Path_new = dataFilePath_new+"Seperated Structures/";
% Path_hybrid = dataFilePath_hybrid+"Seperated Structures/";
% trackTablePath = dataFilePath+dir(fullfile(dataFilePath,"*.trakTable")).name;
% trackTablePath_new = dataFilePath_new+dir(fullfile(dataFilePath_new,"*.trakTable")).name;
% trackTablePath_hybrid = dataFilePath_hybrid+dir(fullfile(dataFilePath_hybrid,"*.trakTable")).name;
% % The path of source file (.nc file)
% srcData="../../MexicoBay/src/MexicoBayData.nc4";
% owDataFilePath="../../MexicoBay/OW/FT_result/";
% Dataset = "MexicoBay";
% 
% 
% property.x = "Longitude";
% property.y = "Latitude";
% property.z = "Depth";
% property.u = "u";
% property.v = "v";
% property.eta = "";
% property.temp = "temperature";
% property.salt = "salinity";
% property.time = "MT";
% 
% 
% ncid = netcdf.open(srcData);
% 
% x_val = ncread(srcData, property.x);
% y_val = ncread(srcData, property.y);
% z_val = ncread(srcData, property.z);
% time_val = ncread(srcData, property.time);
% 
% load("MexicoBay_eddy_clk_data_hybrid.mat", "clockEddy_hybrid");
% load("MexicoBay_eddy_conclk_data_hybrid.mat", "conterclockEddy_hybrid");
% load("MexicoBay_eddy_data_hybrid.mat", "allEddy_hybrid");
% load("MexicoBay_eddy_data_oldWindingAngle.mat", "allEddy");
% load("MexicoBay_eddy_data_newWindingAngle.mat", "allEddy_new");
% 
% 
% % load("MexicoBay_eddy_clk_data.mat");
% % load("MexicoBay_eddy_conclk_data.mat");
% % load("MexicoBay_eddy_data.mat");
% % load("MexicoBay_eddy_path.mat");
% % load("MexicoBay_eddy_history.mat");
% % load("MexicoBay_eddy_graph.mat");
% 
% % save("MexicoBay_eddy_clk_data_hybrid.mat", "clockEddy_hybrid");
% % save("MexicoBay_eddy_conclk_data_hybrid.mat", "conterclockEddy_hybrid");
% % save("MexicoBay_eddy_data_hybrid.mat", "allEddy_hybrid");
% % save("MexicoBay_eddy_path_hybrid.mat", "eddyPath_hybrid");
% % save("MexicoBay_eddy_history_hybrid.mat", "eddyHistory_hybrid");
% % save("MexicoBay_eddy_graph_hybrid.mat", "G_hybrid");
% % 
% % 
% % save("MexicoBay_eddy_data_oldWindingAngle.mat", "allEddy");
% % save("MexicoBay_eddy_path_oldWindingAngle.mat", "eddyPath");
% % save("MexicoBay_eddy_history_oldWindingAngle.mat", "eddyHistory");
% % save("MexicoBay_eddy_graph_oldWindingAngle.mat", "G");
% % 
% % 
% % save("MexicoBay_eddy_data_newWindingAngle.mat", "allEddy_new");
% % save("MexicoBay_eddy_path_newWindingAngle.mat", "eddyPath_new");
% % save("MexicoBay_eddy_history_newWindingAngle.mat", "eddyHistory_new");
% % save("MexicoBay_eddy_graph_newWindingAngle.mat", "G_new");

%% SoutheastAtlantic
% dataFilePath = "../../SEA/windingAngle/FT_result_oldWindingAngle/"; % Output from the windingAngle Eddy Detection
% dataFilePath_new = "../../SEA/windingAngle/FT_result_newWindingAngle/"; % Output from the windingAngle Eddy Detection
% dataFilePath_hybrid = "../../SEA/hybrid/FT_result/";
% 
% Path = dataFilePath+"Seperated Structures/";
% Path_new = dataFilePath_new+"Seperated Structures/";
% Path_hybrid = dataFilePath_hybrid+"Seperated Structures/";
% trackTablePath = dataFilePath+dir(fullfile(dataFilePath,"*.trakTable")).name;
% trackTablePath_new = dataFilePath_new+dir(fullfile(dataFilePath_new,"*.trakTable")).name;
% trackTablePath_hybrid = dataFilePath_hybrid+dir(fullfile(dataFilePath_hybrid,"*.trakTable")).name;
% % The path of source file (.nc file)
% srcData="../../SEA/src/19970101.ocean_5day.nc";
% owDataFilePath="../../RedSea/OW/FT_result/1/";
% Dataset = "SEA";
% 
% 
% property.x = "xh";
% property.y = "yh";
% property.z = "z_l";
% property.u = "u";
% property.v = "v";
% property.eta = "ssh";
% property.temp = "temp";
% property.salt = "salt";
% property.time = "time";
% 
% 
% ncid = netcdf.open(srcData);
% 
% x_val = ncread(srcData, property.x);
% y_val = ncread(srcData, property.y);
% z_val = ncread(srcData, property.z);
% time_val = ncread(srcData, property.time);
% 
% u_val_test = ncread(srcData, property.u, [1,1,1,1], [length(x_val), length(y_val),1,1]);
% v_val_test = ncread(srcData, property.v, [1,1,1,1], [length(x_val), length(y_val),1,1]);
% 
% 
% load("SEA_eddy_clk_data_hybrid.mat", "clockEddy_hybrid");
% load("SEA_eddy_conclk_data_hybrid.mat", "conterclockEddy_hybrid");
% load("SEA_eddy_data_hybrid.mat", "allEddy_hybrid");
% load("SEA_eddy_data_oldWindingAngle.mat", "allEddy");
% load("SEA_eddy_data_newWindingAngle.mat", "allEddy_new");
% 
% 
% % load("SEA_eddy_clk_data.mat");
% % load("SEA_eddy_conclk_data.mat");
% % load("SEA_eddy_data.mat");
% % load("SEA_eddy_path.mat");
% % load("SEA_eddy_history.mat");
% % load("SEA_eddy_graph.mat");
% % 
% % save("SEA_eddy_clk_data_hybrid.mat", "clockEddy_hybrid");
% % save("SEA_eddy_conclk_data_hybrid.mat", "conterclockEddy_hybrid");
% % save("SEA_eddy_data_hybrid.mat", "allEddy_hybrid");
% % save("SEA_eddy_path_hybrid.mat", "eddyPath_hybrid");
% % save("SEA_eddy_history_hybrid.mat", "eddyHistory_hybrid");
% % save("SEA_eddy_graph_hybrid.mat", "G_hybrid");
% % 
% % 
% % save("SEA_eddy_data_oldWindingAngle.mat", "allEddy");
% % save("SEA_eddy_path_oldWindingAngle.mat", "eddyPath");
% % save("SEA_eddy_history_oldWindingAngle.mat", "eddyHistory");
% % save("SEA_eddy_graph_oldWindingAngle.mat", "G");
% % 
% % 
% % save("SEA_eddy_data_newWindingAngle.mat", "allEddy_new");
% % save("SEA_eddy_path_newWindingAngle.mat", "eddyPath_new");
% % save("SEA_eddy_history_newWindingAngle.mat", "eddyHistory_new");
% % save("SEA_eddy_graph_newWindingAngle.mat", "G_new");

%% NorthAtlantic
% dataFilePath = "../../NA/windingAngle/FT_result_oldWindingAngle/"; % Output from the windingAngle Eddy Detection
% dataFilePath_new = "../../NA/windingAngle/FT_result_newWindingAngle/"; % Output from the windingAngle Eddy Detection
% dataFilePath_hybrid = "../../NA/hybrid/FT_result/";
% 
% Path = dataFilePath+"Seperated Structures/";
% Path_new = dataFilePath_new+"Seperated Structures/";
% Path_hybrid = dataFilePath_hybrid+"Seperated Structures/";
% trackTablePath = dataFilePath+dir(fullfile(dataFilePath,"*.trakTable")).name;
% trackTablePath_new = dataFilePath_new+dir(fullfile(dataFilePath_new,"*.trakTable")).name;
% trackTablePath_hybrid = dataFilePath_hybrid+dir(fullfile(dataFilePath_hybrid,"*.trakTable")).name;
% % The path of source file (.nc file)
% srcData="../../NA/src/19960101.ocean_5day.nc";
% owDataFilePath="../../RedSea/OW/FT_result/1/";
% 
% Dataset = "NA";
% 
% property.x = "xh";
% property.y = "yh";
% property.z = "z_l";
% property.u = "u";
% property.v = "v";
% property.eta = "ssh";
% property.temp = "temp";
% property.salt = "salt";
% property.time = "time";
% 
% 
% ncid = netcdf.open(srcData);
% 
% x_val = ncread(srcData, property.x);
% y_val = ncread(srcData, property.y);
% z_val = ncread(srcData, property.z);
% time_val = ncread(srcData, property.time);
% 
% load("NA_eddy_clk_data_hybrid.mat", "clockEddy_hybrid");
% load("NA_eddy_conclk_data_hybrid.mat", "conterclockEddy_hybrid");
% load("NA_eddy_data_hybrid.mat", "allEddy_hybrid");
% load("NA_eddy_data_oldWindingAngle.mat", "allEddy");
% load("NA_eddy_data_newWindingAngle.mat", "allEddy_new");
% 
% % load("NA_eddy_clk_data.mat");
% % load("NA_eddy_conclk_data.mat");
% % load("NA_eddy_data.mat");
% % load("NA_eddy_path.mat");
% % load("NA_eddy_history.mat");
% % load("NA_eddy_graph.mat");
% % 
% % save("NA_eddy_clk_data_hybrid.mat", "clockEddy_hybrid");
% % save("NA_eddy_conclk_data_hybrid.mat", "conterclockEddy_hybrid");
% % save("NA_eddy_data_hybrid.mat", "allEddy_hybrid");
% % save("NA_eddy_path_hybrid.mat", "eddyPath_hybrid");
% % save("NA_eddy_history_hybrid.mat", "eddyHistory_hybrid");
% % save("NA_eddy_graph_hybrid.mat", "G_hybrid");
% % 
% % 
% % save("NA_eddy_data_oldWindingAngle.mat", "allEddy");
% % save("NA_eddy_path_oldWindingAngle.mat", "eddyPath");
% % save("NA_eddy_history_oldWindingAngle.mat", "eddyHistory");
% % save("NA_eddy_graph_oldWindingAngle.mat", "G");
% % 
% % 
% % save("NA_eddy_data_newWindingAngle.mat", "allEddy_new");
% % save("NA_eddy_path_newWindingAngle.mat", "eddyPath_new");
% % save("NA_eddy_history_newWindingAngle.mat", "eddyHistory_new");
% % save("NA_eddy_graph_newWindingAngle.mat", "G_new");


%% red sea data - winding angle
% 
dataFilePath = "../../RedSea/windingAngle/FT_result_oldWindingAngle/1/"; % Output from the windingAngle Eddy Detection
dataFilePath_new = "../../RedSea/windingAngle/FT_result_newWindingAngle/1/"; % Output from the windingAngle Eddy Detection
dataFilePath_hybrid = "../../RedSea/hybrid/FT_result/1/";

Path = dataFilePath+"Seperated Structures/";
Path_new = dataFilePath_new+"Seperated Structures/";
Path_hybrid = dataFilePath_hybrid+"Seperated Structures/";
trackTablePath = dataFilePath+dir(fullfile(dataFilePath,"*.trakTable")).name;
trackTablePath_new = dataFilePath_new+dir(fullfile(dataFilePath_new,"*.trakTable")).name;
trackTablePath_hybrid = dataFilePath_hybrid+dir(fullfile(dataFilePath_hybrid,"*.trakTable")).name;
% The path of source file (.nc file)
srcData="../../RedSea/src/0001/COMBINED_2011013100.nc";
owDataFilePath="../../RedSea/OW/FT_result/1/";
Dataset = "RedSea";


property.x = "XC";
property.y = "YC";
property.z = "Z_MIT40";
property.u = "U";
property.v = "V";
property.eta = "ETA";
property.temp = "TEMP";
property.salt = "SALT";
property.time = "T_AX";


ncid = netcdf.open(srcData);

x_val = ncread(srcData, property.x);
y_val = ncread(srcData, property.y);
z_val = ncread(srcData, property.z);
time_val = ncread(srcData, property.time);

load("Redsea_eddy_clk_data_hybrid.mat", "clockEddy_hybrid");
load("Redsea_eddy_conclk_data_hybrid.mat", "conterclockEddy_hybrid");
load("Redsea_eddy_data_hybrid.mat", "allEddy_hybrid");
load("Redsea_eddy_data_oldWindingAngle.mat", "allEddy");
load("Redsea_eddy_data_newWindingAngle.mat", "allEddy_new");

% load("Redsea_eddy_clk_data.mat");
% load("Redsea_eddy_conclk_data.mat");
% load("Redsea_eddy_data.mat");
% load("Redsea_eddy_path.mat");
% load("Redsea_eddy_history.mat");
% load("Redsea_eddy_graph.mat");
% 
% save("Redsea_eddy_clk_data_hybrid.mat", "clockEddy_hybrid");
% save("Redsea_eddy_conclk_data_hybrid.mat", "conterclockEddy_hybrid");
% save("Redsea_eddy_data_hybrid.mat", "allEddy_hybrid");
% save("Redsea_eddy_path_hybrid.mat", "eddyPath_hybrid");
% save("Redsea_eddy_history_hybrid.mat", "eddyHistory_hybrid");
% save("Redsea_eddy_graph_hybrid.mat", "G_hybrid");
% 
% 
% save("Redsea_eddy_data_oldWindingAngle.mat", "allEddy");
% save("Redsea_eddy_path_oldWindingAngle.mat", "eddyPath");
% save("Redsea_eddy_history_oldWindingAngle.mat", "eddyHistory");
% save("Redsea_eddy_graph_oldWindingAngle.mat", "G");
% 
% 
% save("Redsea_eddy_data_newWindingAngle.mat", "allEddy_new");
% save("Redsea_eddy_path_newWindingAngle.mat", "eddyPath_new");
% save("Redsea_eddy_history_newWindingAngle.mat", "eddyHistory_new");
% save("Redsea_eddy_graph_newWindingAngle.mat", "G_new");



%% Parameters of limitation

spaceLimit.x0 = x_val(1);
spaceLimit.y0 = y_val(1);
spaceLimit.x1 = x_val(end);
spaceLimit.y1 = y_val(end);
livingFrameLimit = 0;
durationLimit = 3;
frameLimit.min = 1;
frameLimit.max = length(ncread(srcData, property.time));
radiusLimit.lower = 3;
radiusLimit.upper = inf;

%% Read Eddy detection information

[allEddy] = readEddyEllipseDetection(Path);
[allEddy_new] = readEddyEllipseDetection(Path_new);
[allEddy_hybrid,clockEddy_hybrid, conterclockEddy_hybrid] = readEddyDetection(Path_hybrid);

%% Read eddy track information

[eddyPath,eddyNodes,G] = eddyActivity(trackTablePath);
[eddyPath_new,eddyNodes_new,G_new] = eddyActivity(trackTablePath_new);
[eddyPath_hybrid,eddyNodes_hybrid,G_hybrid] = eddyActivity(trackTablePath_hybrid);
eddyPathIndex=1:1:length(eddyPath);
eddyPathIndex_new=1:1:length(eddyPath_new);
eddyPathIndex_hybrid=1:1:length(eddyPath_hybrid);
%% Gulf stream extraction and filter
% Get the Gulf stream
% Will be used to distinguish either the eddy it above / below the Gulf
% Stream

% streamRidge = gulfStreamExtraction(srcData,property,spaceLimit);
streamRidge = 0;
northFlag=0;

allEddy = ellipse_propertyFilter(allEddy,spaceLimit,streamRidge, x_val, y_val, northFlag);
allEddy_new = ellipse_propertyFilter(allEddy_new,spaceLimit,streamRidge, x_val, y_val, northFlag);
allEddy_hybrid = propertyFilter(allEddy_hybrid,spaceLimit,streamRidge, x_val, y_val, northFlag);





%% Associate tracking information with detection information
% Get the eddy history variable for an individual eddy. 

for i = 1:1:size(eddyPath,1)
    thisPath=eddyPath(i);
    FindAEddy(i)=cellfun(@(y) contains("1.17",y),thisPath);
end
eddyExists=find(FindAEddy==1);


% eddyPathNum=118;
% eddyPathIndex=1:1:length(eddyPath);
% eddyPathNum=11;
% eddyPathNum=17089;
% eddyHistory = eddyPathProcess(clockEddy,conterclockEddy,eddyPath,eddyPathNum,G,1);






eddyHistory = eddyEllipsePathProcess(allEddy,eddyPath, ...
    G,livingFrameLimit);
eddyHistory_new = eddyEllipsePathProcess(allEddy_new,eddyPath_new, ...
    G_new,livingFrameLimit);
eddyHistory_hybrid = eddyPathProcess(clockEddy_hybrid,conterclockEddy_hybrid,eddyPath_hybrid, ...
    G_hybrid,livingFrameLimit);


eddyRepresentative = cell2mat(cellfun(@(x) x(1,:), eddyHistory,'UniformOutput', false));
eddyRepresentative_new = cell2mat(cellfun(@(x) x(1,:), eddyHistory_new,'UniformOutput', false));
eddyRepresentative_hybrid = cell2mat(cellfun(@(x) x(1,:), eddyHistory_hybrid,'UniformOutput', false));
%% Get eddy statistics plots
% Get the individual eddy statistics information. 

EddyNumRecord = eddyStatFunc(clockEddy,conterclockEddy,eddyHistory,srcData, ...
    property,durationLimit, frameLimit, z_val);

%%
% You could use this to visualize get global statistic figures
patchSize = 30;
frameNum=1;
eddyGlobalStat(allEddy,frameLimit,srcData,radiusLimit,patchSize, property,frameNum);

%% Visualize eddy vertical clip
% Show the vertical clip of the data
eddyVis_vertClip(eddyPathNum,G,eddyPath,eddyHistory,dataFilePath,srcData, property);

%% Visualize individual eddy
% You could use this to visualize an individual eddy path
% eddyIndividualVisIndex = 13016;
eddyIndividualVisIndex = 10334;
% eddyIndividualVisIndex = 9802;
% eddyIndividualVisIndex = 7635;
stretchMode = 1;
eddyIndividualVis(eddyPathIndex,G,eddyPath,eddyHistory,dataFilePath,srcData, property, eddyIndividualVisIndex,stretchMode);

%% Visualize Global eddy statistical property
% You could use this to visualize global distribution

eddyGlobalVis(allEddy, 1, srcData, radiusLimit, property);
%% Visualize path
% You could use this to visualize an all eddies
% You need to change the variable name for the NC file inside the function
% This is the function that generate videos

eddyVis(eddyPathIndex,G,eddyPath,eddyHistory,dataFilePath,srcData, property);


% for i = 1:1:length(EddyNumRecord)
%     eddyPathNum = EddyNumRecord(i);
%     eddyHistory = eddyPathProcess(clockEddy,conterclockEddy,eddyPath, eddyPathNum,G,sizeLimit);
%     eddyVis_2D(eddyPathNum,G,eddyPath,eddyHistory,dataFilePath,srcData, property);
% end


% eddyVis_2D(eddyPathIndex,G,eddyPath,eddyHistory,dataFilePath,srcData, property);

%% Eddy general global vis



%% Eddy global comparison - horizontal
frameIndex = 1:1:length(time_val);
% frameIndex = 1;
depthIndex = 1;

eddySummary_oldWindingAngle = load("/home/weiping/StorageDisk/"+Dataset+"/windingAngle/OriginalData/"+Dataset+"_originalWindingAngle_eddySummary.mat").eddySummary_originalWindingAngle;
eddySummary_newWindingAngle = load("/home/weiping/StorageDisk/"+Dataset+"/windingAngle/OriginalData/"+Dataset+"_newWindingAngle_eddySummary.mat").eddySummary_newWindingAngle;














[distanceSum, distanceRatio,areaSum,areaRatio, statMatrix] = eddyHorizontalMethodComparison_Stat(owDataFilePath,srcData, ...
    eddySummary_oldWindingAngle, eddySummary_newWindingAngle,allEddy_hybrid, ...
    frameIndex, property,depthIndex);

drawDifferenceMap(statMatrix, srcData, property);

%% Eddy temporal comparison

stretchMode = 2;
eddyNumber=15;
eddyNumber_new = 16;
eddyNumber_hybrid = 10;
eddyTemporalComparison(srcData, property, eddyNumber, eddyHistory,dataFilePath, stretchMode);
eddyTemporalComparison(srcData, property, eddyNumber_new, eddyHistory_new,dataFilePath_new, stretchMode);
eddyTemporalComparison(srcData, property, eddyNumber_hybrid, eddyHistory_hybrid,dataFilePath_hybrid, stretchMode);
%% Create sea bed for unity Height map
% unityHeightMap(srcData);

centroid_test = allEddy(:,[1,2,1,2,15,16]);
centroid_test(:,1) = (centroid_test(:,1)-30)/0.04;
centroid_test(:,2) = (centroid_test(:,2)-10)/0.04;
 
%% Eddy interpolation

close all;


interpolationRecord = [];
% for eddyNumber=41
% for eddyNumber=1:1:length(eddyHistory)
for eddyNumber=13
% for eddyNumber=10063

    interpResult = eddyEllipseInterpolation(eddyHistory,dataFilePath,srcData, property, eddyNumber);
    interpolationRecord = [interpolationRecord;interpResult];


    
end
%%
close all
interpolationRecord = load("windingAngleInterpolationRecord.mat").interpolationRecord;

stretchMode = 2;
eddyNumber=13;
stride4Test=2;
ellipseInterpolationVis(srcData, property, eddyNumber, eddyHistory,dataFilePath, interpolationRecord, stride4Test, stretchMode);


interpolationEvalutation = cell2table(interpolationRecord,...
    'VariableNames',{'EddyHistoryIndex','total stride','start frame', 'interpolated frame', 'end frame', 'depth in ground truth', ...
    'Euclidean distance for radius', 'cosine distance (after scaled) for radius', 'correlation for radius', 'SSIM for radius', 'DTW distance','startFrame radius', 'endFrame radius', 'radius from ground truth', 'radius from interpolation', 'groundtruth X','groundtruth Y','interpolated X','interpolated Y','Radius outlier value(abs distance from the mean radius)'});


save("./interpolationEvaluation.mat","interpolationEvalutation");

load("./interpolationEvaluation.mat");


%% long-short eddy sudden change
% windingAngleSummary = load("E:\work\Rutgers\eddy interpolation\RedSea_windingAngle_eddySummary_2frames.mat").eddySummary;
% %% cluster winding angle 3D points
% 
% for timeIndex = 1:1:1
% % for timeIndex = 1:1:length(windingAngleSummary)
%     windingAngleEddySummrayInFrame = windingAngleSummary{timeIndex};
%     windingAngleTracks = cluster3Dpoints(windingAngleEddySummrayInFrame, 0.1,z_val);
% end

%% visualize winding angle clustering
figure,

for trackIndex = 1:length(windingAngleTracks)
    currentTrack = windingAngleTracksamarelData{trackIndex};
    plot3(currentTrack(1,:), currentTrack(2,:), currentTrack(3,:));
    hold on
end
grid on;
set(gca,"ZDir", "reverse");
title("Center of eddies over layers from winding angle approach");

%% match winding angle with hybrid on vertical layers
hybridTracks = allEddy(allEddy(:,15)==1,:);
matches = matchTracks(windingAngleTracks, hybridTracks, 1);




%% long-short eddy sudden change

% 1 - centerX
% 2 - centerY
% 3 - depth (layer)
% 4 - eccentricity
% 5 - longAxis
% 6 - shortAxis






%% comparing between winding angle and hyrbid eddy detection
v = VideoWriter('Hybrid_WindingAngle_Comparison.mp4');
v.FrameRate=1;
v.Quality=100;
windingAngleSummary = load("/home/weiping/StorageDisk/RedSea/windingAngle/amarelData/RedSea_newWindingAngle_eddySummary.mat").eddySummary_newWindingAngle;
saveDir = "/home/weiping/StorageDisk/RedSea/windingAngle/data_afterClustering_Amarel_new/";
dataResolution = 0.04;
track_depthThres = 5;
track_horizontalDistThres = 0.1;
open(v);
for timeIndex = 1
% for timeIndex = 1:1:length(time_val)
% for timeIndex = 1:1:length(windingAngleSummary)
    windingAngleEddySummrayInFrame = windingAngleSummary{timeIndex};
    [windingAngleTracks, eddy_ellipse_objects] = cluster3Dpoints(windingAngleEddySummrayInFrame, track_horizontalDistThres,track_depthThres,x_val,y_val,z_val,dataResolution);

    % match winding angle with hybrid on vertical layers
    hybridTracks = allEddy(allEddy(:,15)==timeIndex,:);
    [matches,unMatchWindingAngle,unMatchHybrid] = matchTracks(windingAngleTracks, hybridTracks, 1.5);
    % winding angle - hyrbid match visualization
    
    VisWindingAngleHybrid(srcData, property, dataFilePath, allEddy, windingAngleSummary, windingAngleTracks, matches, unMatchHybrid,timeIndex,eddy_ellipse_objects);
    F1=getframe(gcf);
    writeVideo(v,F1);
end
close(v);
%% winding angle visualization - 3D comparison

% VisWindingAngleHybrid(srcData, property, dataFilePath, allEddy, windingAngleSummary, matches);

% eddySuddenChangeStat(srcData, property, dataFilePath, allEddy, windingAngleSummary);


eddyPathIndex = 1:1:size(eddyHistory);%(old ones)

eddyPathIndex_new = 1:1:size(eddyHistory_new);%(new ones)

% eddyWindingAngleVis(srcData, property, dataFilePath, allEddy,eddyHistory, ...
%     eddyPathIndex);
eddyWindingAngleVis_twoComparsion(srcData, property, dataFilePath, dataFilePath_new,...
    allEddy,allEddy_new,eddyHistory, eddyHistory_new, eddyPathIndex, eddyPathIndex_new);

xlim([10.5,13]);
ylim([43.5, 45.4]);
zlim([1,1000]);
%% winding angle - hyrbid match visualization

% VisWindingAngleHybrid(srcData, property, dataFilePath, allEddy, windingAngleSummary, matches);

% eddySuddenChangeStat(srcData, property, dataFilePath, allEddy, windingAngleSummary);


eddyPathIndex = 1;

eddySuddenChangeVis(srcData, property, dataFilePath, allEddy,eddyHistory, ...
    eddyPathIndex);
%% Visualize Velocity field at certain layer

timestep = 7;
z_level = 14;
xRange = [200,500];
yRange = [200,500];

u_test = ncread(srcData, property.u, [1,1,z_level,7], [500,500,1,1]);
v_test = ncread(srcData, property.v, [1,1,z_level,7], [500,500,1,1]);

quiver(x_val, y_val, u_test',v_test', 'AutoScaleFactor',3);
daspect([1 1 1])

%% Eddy time series visualization


clkRows = find(clkEddyStat.movingDistance>500);
clkEddyID = table2array(clkEddyStat(clkRows,"id"));

conclkRows = find(conclkEddyStat.movingDistance>500);
conclkEddyID = table2array(conclkEddyStat(conclkRows,"id"));

eddyIndex = [3670,5069,10024,12433];
eddyIndex = [clkEddyID; conclkEddyID];
eddyIndex = [13016];
eddyIndex = [5069];
eddyTimeSeries3DVis(srcData, property, eddyIndex, eddyHistory,dataFilePath);

%% linearProfile
% x_lower = -13.7;
% x_upper = -11.5;
% y_lower = -27.3;
% y_upper = -25.3;

x_lower = 43.5;
x_upper = 46;
y_lower = 10.25;
y_upper = 12.75;

FrameNum = 1;
plotHybridOnly_flag = 0;
fixed_linearProfile_flag = 0;

eddySummary_oldWindingAngle = load("/home/weiping/StorageDisk/"+Dataset+"/windingAngle/OriginalData/"+Dataset+"_originalWindingAngle_eddySummary.mat").eddySummary_originalWindingAngle;
eddySummary_newWindingAngle = load("/home/weiping/StorageDisk/"+Dataset+"/windingAngle/OriginalData/"+Dataset+"_newWindingAngle_eddySummary.mat").eddySummary_newWindingAngle;
linearProfile(srcData, property, x_lower, x_upper, y_lower, y_upper, ...
    FrameNum,allEddy_hybrid,eddySummary_oldWindingAngle, eddySummary_newWindingAngle,plotHybridOnly_flag,fixed_linearProfile_flag)


 %% horizontal attr vis
% % x_lower = -13.7;
% % x_upper = -11.5;
% % y_lower = -27.3;
% % y_upper = -25.3;
% 
% x_lower = 43.5;
% x_upper = 46;
% y_lower = 10.25;
% y_upper = 12.75;
% 
% u_val = ncread(srcData, property.u, [1,1,1,1],[length(x_val),length(y_val),length(z_val),1]);
% v_val = ncread(srcData, property.v, [1,1,1,1],[length(x_val),length(y_val),length(z_val),1]);
% ow_val_prev = calcOW(length(x_val), length(y_val), u_val, v_val);
% 
% figure,
% quiver(x_val(1:3:end), y_val(1:3:end),u_val(1:3:end,1:3:end,1)',v_val(1:3:end,1:3:end,1)',3,"DisplayName","velocity vectors");
% 
% 
% daspect([1 1 1])
