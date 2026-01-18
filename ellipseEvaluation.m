function [results] = ellipseEvaluation(centerX_gt, centerY_gt,lengthX_gt, lengthY_gt, theta_gt, centerX_interp, centerY_interp,lengthX_interp, lengthY_interp, theta_interp)
%ELLIPSEEVALUATION Summary of this function goes here
%   Detailed explanation goes here


% centerX_evaluation = eddyCoordEvaluationFunc(centerX_gd, centerX_interp(2,:)');
% centerY_evaluation = eddyCoordEvaluationFunc(centerY_gd, centerY_interp(2,:)');
% theta_evaluation = eddyAngleEvaluationFunc(theta_gd', theta_interp(2,:));
% lengthX_evaluation = eddyScalarEvaluationFunc(lengthX_gd, lengthX_interp(2,:)');
% lengthY_evaluation = eddyScalarEvaluationFunc(lengthY_gd, lengthY_interp(2,:)');

depth_gt = 1:1:length(theta_gt);
depth_interp = 1:1:length(theta_interp);
pca_components = 3;

% ---- Step 1: 角度转换为向量形式 ----
angle_gt = [cosd(theta_gt(:))'; sind(theta_gt(:))'];
angle_interp = [cosd(theta_interp(:))'; sind(theta_interp(:))'];

% ---- Step 2: 统一深度范围并下采样到共同深度 ----
min_len = min(length(depth_gt), length(depth_interp));
depth_common = linspace(max(min(depth_gt), min(depth_interp)), ...
                        min(max(depth_gt), max(depth_interp)), min_len);


data_gt = [
    centerX_gt(:)';
    centerY_gt(:)';
    angle_gt(1,:);
    angle_gt(2,:);
    lengthX_gt(:)';
    lengthY_gt(:)'
];

data_interp = [
    centerX_interp(:)';
    centerY_interp(:)';
    angle_interp(1,:);
    angle_interp(2,:);
    lengthX_interp(:)';
    lengthY_interp(:)'
];



% 安全插值函数（带外推）
interp1_safe = @(x, y) interp1(x, y, depth_common, 'linear', 'extrap');

% 构造插值后的变量矩阵（每行一个变量）
vars_gt = [
    interp1_safe(depth_gt, centerX_gt(:)');
    interp1_safe(depth_gt, centerY_gt(:)');
    interp1_safe(depth_gt, angle_gt(1,:));
    interp1_safe(depth_gt, angle_gt(2,:));
    interp1_safe(depth_gt, lengthX_gt(:)');
    interp1_safe(depth_gt, lengthY_gt(:)')
];

vars_interp = [
    interp1_safe(depth_interp, centerX_interp(:)');
    interp1_safe(depth_interp, centerY_interp(:)');
    interp1_safe(depth_interp, angle_interp(1,:));
    interp1_safe(depth_interp, angle_interp(2,:));
    interp1_safe(depth_interp, lengthX_interp(:)');
    interp1_safe(depth_interp, lengthY_interp(:)')
];

% ---- Step 3: PCA + RMSE + Cosine Similarity ----
pca_components = 3;  % 可调参数，主成分数量
X_all = [vars_gt'; vars_interp'];
X_all = (X_all - mean(X_all)) ./ std(X_all);
[~, score] = pca(X_all);
score_gt = score(1:min_len, 1:pca_components);
score_interp = score(min_len+1:end, 1:pca_components);

results = struct();
results.PCA_RMSE = sqrt(mean((score_gt - score_interp).^2, 'all'));
results.PCA_CosSim = dot(score_gt(:), score_interp(:)) / (norm(score_gt(:)) * norm(score_interp(:)));


labels = {'centerX', 'centerY', 'angleX', 'angleY', 'lengthX', 'lengthY'};

% ---- Step 5: DTW ----
results.DTW_Dist = zeros(1, length(labels));
for i = 1:size(data_gt, 1)
    results.DTW_Dist(i) = dtw(data_gt(i,:)', data_interp(i,:)');
end

results.labels = labels;



end

