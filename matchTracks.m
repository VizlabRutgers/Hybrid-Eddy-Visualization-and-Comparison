function [matches,unMatchWindingAngle,unMatchHybrid] = matchTracks(windingAngleTracks, hybridTracks, maxDist)
% tracksA, tracksB: cell arrays of size M, N；每个元素为 n×3 矩阵
% maxDist: 最大允许匹配距离阈值
% matches: K×2 数组，每行 [i, j] 表示 A(i) 匹配 B(j)

if nargin < 3
    maxDist = 5;  % 根据场景调整
end

M = numel(windingAngleTracks);
N = size(hybridTracks,1);

% 1. 计算质心

centroidsWindingAngle = cellfun(@(x) x(1:2,1), windingAngleTracks, 'UniformOutput', false);
centroidsWindingAngle = cellfun(@(x) cell2mat(x), centroidsWindingAngle, 'UniformOutput', false);
centroidsWindingAngle = cell2mat(centroidsWindingAngle)';
centroidsHybrid = hybridTracks(:,1:2);


% 2. 构造距离矩阵
D = pdist2(centroidsWindingAngle, centroidsHybrid);  % M×N

% 4. 匈牙利算法求解最小分配
% MATLAB 自带：matchpairs （需要 Statistics and Machine Learning Toolbox）
[matches,unMatchWindingAngle,unMatchHybrid] = matchpairs(D, maxDist);


end
