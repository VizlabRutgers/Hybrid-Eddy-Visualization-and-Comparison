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

%% na data
% dataFilePath = "/home/weiping/StorageDisk/NA_data_highResolution/hybrid/FT_result/"; % Output from the Hybrid Eddy Detection
% Path = dataFilePath+"Seperated Structures/";
% trackTablePath = dataFilePath+dir(fullfile(dataFilePath,"*.trakTable")).name;
% % The path of source file (.nc file)
% 
% % To open a remote dataset, use its URL:
% srcData='/home/weiping/StorageDisk/NA_data_highResolution/src/19960101.ocean_5day.nc';
% owDataFilePath="../../FT_OW/1/";
% 
% 
% 
% property.x = "xh";
% property.y = "yh";
% property.z = "z_l";
% property.u = "u";
% property.v = "v";
% property.eta = "nan";
% property.temp = "temp";
% property.salt = "salt";
% property.time = "time";
% 
% ncid = netcdf.open (srcData);
% 
% x_val = ncread(srcData, property.x);
% y_val = ncread(srcData, property.y);
% z_val = ncread(srcData, property.z);
% 
% load("NA_eddy_clk_data.mat");
% load("NA_eddy_conclk_data.mat");
% load("NA_eddy_data");
% load("NA_eddy_path.mat");
% load("NA_eddy_history.mat");
% load("NA_eddy_graph.mat");
% % eddyPathIndex=1:1:length(eddyPath);
% 
% % save("NA_eddy_clk_data.mat", "clockEddy");
% % save("NA_eddy_conclk_data.mat", "conterclockEddy");
% % save("NA_eddy_data.mat", "allEddy");
% % save("NA_eddy_path.mat", "eddyPath");
% % save("NA_eddy_history.mat", "eddyHistory");
% % save("NA_eddy_graph.mat", "G");
% % eddyPathIndex=1:1:length(eddyPath);
% 
% % 
% % 
% % eddyStatistic_NA = array2table(allEddy(:,1:16),...
% %     'VariableNames',{'eddyCenter_x (degree)','eddyCenter_y (degree)','first_point_on_boundary_x', 'first_point_on_boundary_y', 'first_point_on_boundary_z(meter)', ...
% %     'OW (first point)', 'velocity_u (first point)', 'velocity_v (first point)', 'vorticity (first point)', 'salinity (first point)', 'temperature (first point)', 'all zeros (leave for future use)','radius on surface(pixel)', 'clockwiseFlage (clockwise = 1)', 'frameIndex', 'eddyIndex in current frame'});
% % 


%% red sea data
% % dataFilePath = "/home/weiping/data/ft_changes/FT_test/1/"; % Output from the Hybrid Eddy Detection
% dataFilePath = "/home/weiping/StorageDisk/RedSea/hybrid/FT_result/1/"; % Output from the Hybrid Eddy Detection
% Path = dataFilePath+"Seperated Structures/";
% trackTablePath = dataFilePath+dir(fullfile(dataFilePath,"*.trakTable")).name;
% % The path of source file (.nc file)
% srcData="/home/weiping/data/SciViz/SciVisContest2020/ensembles/0001/COMBINED_2011013100.nc";
% owDataFilePath="../../FT_OW/1/";
% 
% property.x = "XC";
% property.y = "YC";
% property.z = "Z_MIT40";
% property.u = "U";
% property.v = "V";
% property.eta = "ETA";
% property.temp = "TEMP";
% property.salt = "SALT";
% property.time = "T_AX";
% 
% 
% ncid = netcdf.open(srcData);
% 
% x_val = ncread(srcData, property.x);
% y_val = ncread(srcData, property.y);
% z_val = ncread(srcData, property.z);
% time_val = ncread(srcData, property.time);
% 
% load("Redsea_eddy_clk_data.mat");
% load("Redsea_eddy_conclk_data.mat");
% load("Redsea_eddy_data.mat");
% load("Redsea_eddy_path.mat");
% load("Redsea_eddy_history.mat");
% load("Redsea_eddy_graph.mat");
% eddyPathIndex=1:1:length(eddyPath);
% % 
% % save("Redsea_eddy_clk_data.mat", "clockEddy");
% % save("Redsea_eddy_conclk_data.mat", "conterclockEddy");
% % save("Redsea_eddy_data.mat", "allEddy");
% % save("Redsea_eddy_path.mat", "eddyPath");
% % save("Redsea_eddy_history.mat", "eddyHistory");
% % save("Redsea_eddy_graph.mat", "G");
% % eddyPathIndex=1:1:length(eddyPath);

%% Mexico Bay

% dataFilePath = "/home/weiping/data/ft_changes/FT_test/1/"; % Output from the Hybrid Eddy Detection
dataFilePath = "/home/weiping/StorageDisk/Mexico_Bay/hybrid/FT_result/"; % Output from the Hybrid Eddy Detection
Path = dataFilePath+"Seperated Structures/";
trackTablePath = dataFilePath+dir(fullfile(dataFilePath,"*.trakTable")).name;
% The path of source file (.nc file)
srcData="/home/weiping/StorageDisk/Mexico_Bay/MexicoBayData.nc4";
% owDataFilePath="../../FT_OW/1/";

property.x = "Longitude";
property.y = "Latitude";
property.z = "Depth";
property.u = "u";
property.v = "v";
property.eta = "nan";
property.temp = "temperature";
property.salt = "salinity";
property.time = "MT";


ncid = netcdf.open(srcData);

x_val = ncread(srcData, property.x);
y_val = ncread(srcData, property.y);
z_val = ncread(srcData, property.z);
time_val = ncread(srcData, property.time);

% load("MexicoBay_eddy_clk_data.mat");
% load("MexicoBay_eddy_conclk_data.mat");
% load("MexicoBay_eddy_data.mat");
% load("MexicoBay_eddy_path.mat");
% load("MexicoBay_eddy_history.mat");
% load("MexicoBay_eddy_graph.mat");
% eddyPathIndex=1:1:length(eddyPath);
% 
% save("Redsea_eddy_clk_data.mat", "clockEddy");
% save("Redsea_eddy_conclk_data.mat", "conterclockEddy");
% save("Redsea_eddy_data.mat", "allEddy");
% save("Redsea_eddy_path.mat", "eddyPath");
% save("Redsea_eddy_history.mat", "eddyHistory");
% save("Redsea_eddy_graph.mat", "G");
% eddyPathIndex=1:1:length(eddyPath);


%% Parameters of limitation

spaceLimit.x0 = x_val(1);
spaceLimit.y0 = y_val(1);
spaceLimit.x1 = x_val(end);
spaceLimit.y1 = y_val(end);
livingFrameLimit = 0;
durationLimit = 3;
frameLimit.min = 1;
frameLimit.max = inf;
radiusLimit.lower = 3;
radiusLimit.upper = inf;


%% Pre-processing (not necessary if you load data before)
% % Read Eddy detection information
% 
[allEddy, clockEddy, conterclockEddy] = readEddyDetection(Path);

% Read eddy track information

[eddyPath,eddyNodes,G] = eddyActivity(trackTablePath);
eddyPathIndex=1:1:length(eddyPath);
% Gulf stream extraction and filter
% Get the Gulf stream
% Will be used to distinguish either the eddy it above / below the Gulf
% Stream

% streamRidge = gulfStreamExtraction(srcData,property,spaceLimit);
streamRidge = 0;
northFlag=0;

allEddy = propertyFilter(allEddy,spaceLimit,radiusLimit, streamRidge, x_val, y_val, northFlag);
clockEddy = propertyFilter(clockEddy,spaceLimit,radiusLimit, streamRidge, x_val, y_val, northFlag);
conterclockEddy = propertyFilter(conterclockEddy,spaceLimit, radiusLimit, streamRidge, x_val, y_val, northFlag);




% Associate tracking information with detection information
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






eddyHistory = eddyPathProcess(clockEddy,conterclockEddy,eddyPath, ...
    G,livingFrameLimit);




%% Global eddy movement statistics plots
% Get the global eddy path statistics information. 
% Including movement distance, existing time, 

[clkEddyStat, conclkEddyStat, EddyNumRecord ] = eddyStatFunc(clockEddy,conterclockEddy,eddyHistory,srcData, ...
    property,durationLimit, frameLimit, z_val);

%% Eddy path statistics
% Analyze the statistical change of the eddy over the whole path
eddyPathStatResult = eddyPathStat(eddyHistory,srcData, property,dataFilePath);
eddyPathStatResult = eddyPathStatResult';
eddyPathStatWhole = {};
for i = 1:1: length(eddyPathStatResult)
    eddyPathStatWhole = [eddyPathStatWhole; table2cell(eddyPathStatResult{i})];
end

% figure,
% scatter(cell2mat(eddyPathStatWhole(:,5)), cell2mat(eddyPathStatWhole(:,4)));
% xlabel("surface radius (pixel)");
% ylabel("depth (meter)");
% 
% figure,
% scatter(cell2mat(eddyPathStatWhole(:,6)), cell2mat(eddyPathStatWhole(:,4)));
% xlabel("volume");
% ylabel("depth (meter)");
% 
% figure,
% scatter(cell2mat(eddyPathStatWhole(:,7)), cell2mat(eddyPathStatWhole(:,4)));
% xlabel("minimum OW");
% ylabel("depth (meter)");
% 
% figure,
% scatter(cell2mat(eddyPathStatWhole(:,8)), cell2mat(eddyPathStatWhole(:,4)));
% xlabel("maximum angular");
% ylabel("depth (meter)");
% 
% figure,
% scatter(cell2mat(eddyPathStatWhole(:,8)), cell2mat(eddyPathStatWhole(:,7)));
% xlabel("maximum angular");
% ylabel("OW");
%% Eddy global distribution Visualization
% You could use this to visualize global eddy distribution
% Including eddy number, vorticity, size, etc.
% Also use large patches to visualize the average distribution in an area
% if FrameIndex not defined, then use frameLimit
globalVisFrameIndex = 1;
patchSize = 30;
eddyGlobalStat(allEddy,frameLimit,srcData,radiusLimit,patchSize, property, globalVisFrameIndex);

%% Visualize eddy vertical clip
% Show the vertical clip of the data
eddyVis_vertClip(eddyPathNum,G,eddyPath,eddyHistory,dataFilePath,srcData, property);

%% Visualize individual eddy path
% You could use this to visualize an individual eddy path
% eddyIndividualVisIndex = 13016;
% eddyIndividualVisIndex = 10334;
% eddyIndividualVisIndex = 9802;
% eddyIndividualVisIndex = 7635;
% eddyIndividualVisIndex = 9037;
% eddyIndividualVisIndex = 9947;
% eddyIndividualVisIndex = 10150;
% eddyIndividualVisIndex = 2155;
% eddyIndividualVisIndex = 3670;
% eddyIndividualVisIndex = 5069;
eddyIndividualVisIndex = 10024;
eddyIndividualVisIndex = 12433;
stretchMode = 1;


eddyIndividualVisIndex = 11;
eddyIndividualVis(eddyPathIndex,G,eddyPath,eddyHistory,dataFilePath,srcData, property, eddyIndividualVisIndex,stretchMode);

%% Visualize Global eddy
% You could use this to visualize global distribution

eddyGlobalVis(allEddy, 1, srcData, radiusLimit, property,dataFilePath);
%% filter Eddy


testEddy = allEddy(allEddy(:,15) == 73,:);
eddyFrame = 73;
eddyIndex = 376;

eddyHistory_toBeFiltered = cellfun(@(x) x(end,:), eddyHistory);
eddyHistory_toBeFiltered = cell2mat(eddyHistory_toBeFiltered);
[~,loc] = ismember(testEddy(145,:),eddyHistory_toBeFiltered,"rows");

%% Visualize path
% You could use this to visualize an eddy path
% You need to change the variable name for the NC file inside the function
% This is the function that generate videos

eddyPathIndex = 5;
eddyVis(eddyPathIndex,G,eddyPath,eddyHistory,dataFilePath,srcData, property);


% for i = 1:1:length(EddyNumRecord)
%     eddyPathNum = EddyNumRecord(i);
%     eddyHistory = eddyPathProcess(clockEddy,conterclockEddy,eddyPath, eddyPathNum,G,sizeLimit);
%     eddyVis_2D(eddyPathNum,G,eddyPath,eddyHistory,dataFilePath,srcData, property);
% end


% eddyVis_2D(eddyPathIndex,G,eddyPath,eddyHistory,dataFilePath,srcData, property);



%% Eddy Method comparison
% Compare with different methods (Hybrid / OW / winding angle)
frameIndex = 1;
windingAngleDataFilePath = "";
eddyMethodComparison(owDataFilePath,windingAngleDataFilePath,srcData, allEddy, frameIndex, property);


%% Create sea bed for unity Height map
% unityHeightMap(srcData);

centroid_test = allEddy(:,[1,2,1,2,15,16]);
centroid_test(:,1) = (centroid_test(:,1)-30)/0.04;
centroid_test(:,2) = (centroid_test(:,2)-10)/0.04;

%% Self-define
% 
% eddyPathIndex = 11;
% frameIndex = 18:1:19;

eddyPathIndex = 10024;
frameIndex = 36;

selfDefineVis(allEddy, frameIndex, eddyPathIndex, srcData, property,dataFilePath,eddyHistory);

%% Eddy interpolation

close all;

interpolationRecord = [];
% for eddyNumber=41
% for eddyNumber=1:1:length(eddyHistory)
% for eddyNumber=11
for eddyNumber=10024

    evaluationResult = eddyInterpolation(eddyPathIndex,G,eddyPath,eddyHistory,dataFilePath,srcData, property, eddyNumber);
    interpolationRecord = [interpolationRecord; evaluationResult];
    
end

stride4Test = 2;
stretchMode = 2;
interpolationVis(srcData, property, eddyNumber, eddyHistory,dataFilePath, interpolationRecord, stride4Test, stretchMode);





interpolationEvalutation = cell2table(interpolationRecord,...
    'VariableNames',{'EddyHistoryIndex','total stride','start frame', 'interpolated frame', 'end frame', 'depth in ground truth', ...
    'Euclidean distance for radius', 'cosine distance (after scaled) for radius', 'correlation for radius', 'SSIM for radius', 'DTW distance','startFrame radius', 'endFrame radius', 'radius from ground truth', 'radius from interpolation', 'groundtruth X','groundtruth Y','interpolated X','interpolated Y','Radius outlier value(abs distance from the mean radius)'});


save("./interpolationEvaluation.mat","interpolationEvalutation");

load("./interpolationEvaluation.mat");


%% long-short eddy sudden change
% 1 - centerX
% 2 - centerY
% 3 - depth (layer)
% 4 - eccentricity
% 5 - longAxis
% 6 - shortAxis


windingAngleSummary = load("/home/weiping/StorageDisk/Mexico_Bay/windingAngle/OriginalData/MexicoBay_originalWindingAngle_eddySummary.mat").eddySummary_originalWindingAngle;



%% Cluster 3D windingAngle centroids and corresponding volumes
saveDir = "~/StorageDisk/Mexico_Bay/windingAngle/data_afterClustering_originalWindingAngle/";
dataResolution = 0.04;
% for timeIndex = 1:1:2
for timeIndex = 1:1:length(windingAngleSummary)
    windingAngleEddySummrayInFrame = windingAngleSummary{timeIndex};
    [windingAngleTracks, eddy_ellipse_objects] = cluster3Dpoints(windingAngleEddySummrayInFrame, 0.1,x_val,y_val,z_val,dataResolution);
    
    for objectIndex = 1:1:length(eddy_ellipse_objects)
        currentObject = eddy_ellipse_objects{objectIndex};
        writematrix(currentObject, saveDir+"Frame_"+num2str(timeIndex)+"_ellipsoidEddy_"+num2str(objectIndex)+".txt",'Delimiter',' ')  
    end
    
end

%% comparing between winding angle and hyrbid eddy detection
v = VideoWriter('Hybrid_WindingAngle_Comparison.mp4');
v.FrameRate=1;
v.Quality=100;
open(v);
dataResolution = 0.04;
for timeIndex = 1
% for timeIndex = 1:1:length(time_val)
% for timeIndex = 1:1:length(windingAngleSummary)
    windingAngleEddySummrayInFrame = windingAngleSummary{timeIndex};
    [windingAngleTracks, eddy_ellipse_objects] = cluster3Dpoints(windingAngleEddySummrayInFrame, 0.1,x_val,y_val,z_val,dataResolution);

    % match winding angle with hybrid on vertical layers
    hybridTracks = allEddy(allEddy(:,15)==timeIndex,:);
    [matches,unMatchWindingAngle,unMatchHybrid] = matchTracks(windingAngleTracks, hybridTracks, 1.5);
    % winding angle - hyrbid match visualization
    
    VisWindingAngleHybrid(srcData, property, dataFilePath, allEddy, windingAngleSummary, windingAngleTracks, matches, unMatchHybrid,timeIndex,eddy_ellipse_objects);
    F1=getframe(gcf);
    writeVideo(v,F1);
end
close(v);


%% Check eddy sudden change
eddySuddenChangeStat(srcData, property, dataFilePath, allEddy, windingAngleSummary);





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

%% eddy biovis

eddyPathIndex = 10;
eddyBioVis(srcData, property, dataFilePath, allEddy, eddyHistory,eddyPathIndex);

%% eddy biovis surface
timeIndex = 2;

eddyBioVis_surface(srcData, property, dataFilePath, allEddy, eddyHistory,timeIndex);