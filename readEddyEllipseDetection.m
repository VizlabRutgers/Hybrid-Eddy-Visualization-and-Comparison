function [allEddy] = readEddyEllipseDetection(dataFilePath)
%READEDDYDETECTION 此处显示有关此函数的摘要
%   此处显示详细说明

% Definition of each coloumn in the output by sequence

% eddy center coordinates_x on surface
% eddy center coordinates_y on surface
% eddy boundary coordinates_x: this only indicate the first point on the boundary for each eddy
% eddy boundary coordinates_y: this only indicate the first point on the boundary for each eddy
% eddy boundary coordinates_z: this only indicate the first point on the boundary for each eddy



Path = dataFilePath;
file = dir(fullfile(Path,'*.uocd'));
filenames = {file.name};

allEddy=[];
clockEddy = [];
conterclockEddy = [];

% for i=1:1:length(y_val)-1
%     distance_per_pixel(i) = m_lldist([0,0], [y_val(i), y_val(i+1)]);
% end

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
    

    allEddy=[allEddy;[data(1,:),max(data(:,5))]];
end


end

