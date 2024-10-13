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
dataFilePath = "../../NA_dataset/NA_new_1_25/19960101/"; % Output from the Hybrid Eddy Detection
Path = dataFilePath+"Seperated Structures/";
trackTablePath = dataFilePath+dir(fullfile(dataFilePath,"*.trakTable")).name;
% The path of source file (.nc file)

% To open a remote dataset, use its URL:
srcData='http://antares.esm.rutgers.edu:8080/thredds/dodsC/MOM6/ESMG/NWA25/SIM3.0_v2/ocean_5days/19960101.ocean_5day.nc';
owDataFilePath="../../FT_OW/1/";



property.x = "xh";
property.y = "yh";
property.z = "z_l";
property.u = "u";
property.v = "v";
property.eta = "nan";
property.temp = "temp";
property.salt = "salt";
property.time = "time";

ncid = netcdf.open (srcData);

x_val = ncread(srcData, property.x);
y_val = ncread(srcData, property.y);
z_val = ncread(srcData, property.z);

load("NA_eddy_clk_data.mat");
load("NA_eddy_conclk_data.mat");
load("NA_eddy_data");
load("NA_eddy_path.mat");
load("NA_eddy_history.mat");
load("NA_eddy_graph.mat");
eddyPathIndex=1:1:length(eddyPath);

% 
% 
% eddyStatistic_NA = array2table(allEddy(:,1:16),...
%     'VariableNames',{'eddyCenter_x (degree)','eddyCenter_y (degree)','first_point_on_boundary_x', 'first_point_on_boundary_y', 'first_point_on_boundary_z(meter)', ...
%     'OW (first point)', 'velocity_u (first point)', 'velocity_v (first point)', 'vorticity (first point)', 'salinity (first point)', 'temperature (first point)', 'all zeros (leave for future use)','radius on surface(pixel)', 'clockwiseFlage (clockwise = 1)', 'frameIndex', 'eddyIndex in current frame'});
% 


%% red sea data
% dataFilePath = "../../FT_result/1/"; % Output from the Hybrid Eddy Detection
% Path = dataFilePath+"Seperated Structures/";
% trackTablePath = dataFilePath+dir(fullfile(dataFilePath,"*.trakTable")).name;
% % The path of source file (.nc file)
% srcData="../../redsea/0001/COMBINED_2011013100.nc";
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
% 
% load("Redsea_eddy_clk_data.mat");
% load("Redsea_eddy_conclk_data.mat");
% load("Redsea_eddy_data.mat");
% load("Redsea_eddy_path.mat");
% load("Redsea_eddy_history.mat");
% load("Redsea_eddy_graph.mat");
% eddyPathIndex=1:1:length(eddyPath);

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

[allEddy, clockEddy, conterclockEddy] = readEddyDetection(Path);

%% Read eddy track information

[eddyPath,eddyNodes,G] = eddyActivity(trackTablePath);
%% Gulf stream extraction and filter
% Get the Gulf stream
% Will be used to distinguish either the eddy it above / below the Gulf
% Stream

% streamRidge = gulfStreamExtraction(srcData,property,spaceLimit);
streamRidge = 0;
northFlag=0;

allEddy = propertyFilter(allEddy,spaceLimit,radiusLimit, streamRidge, x_val, y_val, northFlag);
clockEddy = propertyFilter(clockEddy,spaceLimit,radiusLimit, streamRidge, x_val, y_val, northFlag);
conterclockEddy = propertyFilter(conterclockEddy,spaceLimit, radiusLimit, streamRidge, x_val, y_val, northFlag);




%% Associate tracking information with detection information
% Get the eddy history variable for an individual eddy. 

for i = 1:1:size(eddyPath,1)
    thisPath=eddyPath(i);
    FindAEddy(i)=cellfun(@(y) contains("1.17",y),thisPath);
end
eddyExists=find(FindAEddy==1);


% eddyPathNum=118;
eddyPathIndex=1:1:length(eddyPath);
% eddyPathNum=11;
% eddyPathNum=17089;
% eddyHistory = eddyPathProcess(clockEddy,conterclockEddy,eddyPath,eddyPathNum,G,1);






eddyHistory = eddyPathProcess(clockEddy,conterclockEddy,eddyPath, ...
    G,livingFrameLimit);


%% Get eddy statistics plots
% Get the individual eddy statistics information. 

EddyNumRecord = eddyStatFunc(clockEddy,conterclockEddy,eddyHistory,srcData, ...
    property,durationLimit, frameLimit, z_val);


% You could use this to visualize get global statistic figures
patchSize = 30;
eddyGlobalStat(allEddy,frameLimit,srcData,radiusLimit,patchSize, property);

%% Visualize eddy vertical clip
% Show the vertical clip of the data
eddyVis_vertClip(eddyPathNum,G,eddyPath,eddyHistory,dataFilePath,srcData, property);

%% Visualize individual eddy
% You could use this to visualize an individual eddy path
eddyIndividualVisIndex = 10063;
eddyIndividualVis(eddyPathIndex,G,eddyPath,eddyHistory,dataFilePath,srcData, property, eddyIndividualVisIndex);

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

%% Eddy global comparison
frameIndex = 1;
windingAngleDataFilePath = "";
eddyMethodComparison(owDataFilePath,windingAngleDataFilePath,srcData, allEddy, frameIndex, property);


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
% for eddyNumber=11
for eddyNumber=10063

    evaluationResult = eddyInterpolation(eddyPathIndex,G,eddyPath,eddyHistory,dataFilePath,srcData, property, eddyNumber);
    interpolationRecord = [interpolationRecord; evaluationResult];


    
end

stride4Test = 2;
stretchMode = 1;
interpolationVis(srcData, property, eddyNumber, eddyHistory,dataFilePath, interpolationRecord, stride4Test, stretchMode);





interpolationEvalutation = cell2table(interpolationRecord,...
    'VariableNames',{'EddyHistoryIndex','total stride','start frame', 'interpolated frame', 'end frame', 'depth in ground truth', ...
    'Euclidean distance for radius', 'cosine distance (after scaled) for radius', 'correlation for radius', 'SSIM for radius', 'DTW distance','startFrame radius', 'endFrame radius', 'radius from ground truth', 'radius from interpolation', 'groundtruth X','groundtruth Y','interpolated X','interpolated Y','Radius outlier value(abs distance from the mean radius)'});


save("./interpolationEvaluation.mat","interpolationEvalutation");

load("./interpolationEvaluation.mat");


