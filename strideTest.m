function [] = strideTest(eddy_forTest, titleStr)
%STRIDETEST Summary of this function goes here
%   Detailed explanation goes here

close all;

eddy_forTest.stride = eddy_forTest.("end frame")-eddy_forTest.("start frame");

figure
boxchart(eddy_forTest.stride, table2array(eddy_forTest(:,"Euclidean distance")));
title("Euclidean distance " + titleStr);
xlabel("stride");

figure,
boxchart(eddy_forTest.stride, table2array(eddy_forTest(:,"cosine distance (after scaled)")));
title("cosine distance (after scaled) "+ titleStr);
xlabel("stride");

figure,
boxchart(eddy_forTest.stride, table2array(eddy_forTest(:,"correlation")));
title("correlation "+ titleStr);
xlabel("stride");

figure,
boxchart(eddy_forTest.stride, table2array(eddy_forTest(:,"SSIM")));
title("SSIM "+ titleStr);
xlabel("stride");

figure,
boxchart(eddy_forTest.stride, table2array(eddy_forTest(:,"DTW distance")));
title("DTW "+ titleStr);
xlabel("stride");

stride = cellstr(num2str(eddy_forTest.stride));
figure
violinplot(table2array(eddy_forTest(:,"Euclidean distance")), stride);
ylabel('Euclidean distance');
xlim([0.5, 7.5]);
ylim([0,1]);

stride = cellstr(num2str(eddy_forTest.stride));
figure
violinplot(table2array(eddy_forTest(:,"DTW distance")), stride);
ylabel('DTW distance');
xlim([0.5, 7.5]);
ylim([0,2]);


end

