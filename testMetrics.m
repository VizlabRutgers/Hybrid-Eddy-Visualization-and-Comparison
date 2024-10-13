function [] = testMetrics(eddy_statistic,titleStr)
%TESTMETRICS Summary of this function goes here
%   Detailed explanation goes here
close all;

figure,
euclidean_box = boxplot(table2array(eddy_statistic(:,7)));
title("Euclidean Distance for " + titleStr);
saveas(gcf,titleStr + '_Euclidean Distance.png')

figure,
cosine_box = boxplot(table2array(eddy_statistic(:,8)));
title("Cosine similarity for " + titleStr);
saveas(gcf,titleStr + '_Cosine similarity.png')

figure,
corr_box = boxplot(table2array(eddy_statistic(:,9)));
title("Correlation for " + titleStr);
saveas(gcf,titleStr + '_Correlation.png')

figure,
ssin_box = boxplot(table2array(eddy_statistic(:,10)));
title("SSIM for " + titleStr);
saveas(gcf,titleStr + '_SSIM.png')

figure,
dtw_box = boxplot(table2array(eddy_statistic(:,15)));
title("dtw for " + titleStr);
saveas(gcf,titleStr + '_dtw.png')

end

