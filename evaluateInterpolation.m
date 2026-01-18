clc;
clear all;


load("./interpolationEvaluation.mat");


eddy_statistic = table2array(interpolationEvalutation);
eddy_13 = eddy_statistic(eddy_statistic(:,1)==13,:);


eddy_13_table = interpolationEvalutation(interpolationEvalutation.EddyIndex == 13,:);


close all;
testMetrics(eddy_13, "eddy 13");
testMetrics(eddy_statistic, "all eddy");

strideTest(eddy_13_table,"eddy 13");
strideTest(interpolationEvalutation, "all eddy");

radius_gd =  table2array(interpolationEvalutation(:,13));
radius_interp =  table2array(interpolationEvalutation(:,14));


stride_Index = table2array(interpolationEvalutation(:,3));

radius_gd_mean = cellfun(@mean,radius_gd);
radius_interp_mean = cellfun(@mean,radius_interp);

stride = cellstr(num2str(stride_Index));
boxplot(abs(radius_gd_mean - radius_interp_mean),stride_Index);


