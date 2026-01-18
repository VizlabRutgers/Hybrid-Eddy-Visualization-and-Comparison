function [isosurface_handle] = owIndividualIsosurface(ax, minX_index, maxX_index, minY_index, maxY_index, dataVolume, voxel_threshold,value_threshold, color, alpha, x_val, y_val, z_val)
%OWINDIVIDUALISOSURFACE 此处显示有关此函数的摘要
%   此处显示详细说明

[gridX, gridY, gridZ] = meshgrid(x_val(minY_index:maxY_index), y_val(minX_index:maxX_index), z_val);
[gridX_2D, gridY_2D] = meshgrid(x_val(minY_index:maxY_index), y_val(minX_index:maxX_index));

data_filtered = dataVolume(minX_index:maxX_index, minY_index:maxY_index,:);
% data_filtered(data_filtered<0) = NaN;
data_filtered = data_filtered<value_threshold;


% Step 2: 连通区域分析
CC = bwconncomp(data_filtered, 26);  % 26 表示3D全连通

% Step 3: 计算每个区域体素数
numPixels = cellfun(@numel, CC.PixelIdxList);

% Step 4: 设置体素数量阈值（如小于100的去掉）
validIdx = find(numPixels >= voxel_threshold);

% Step 5: 创建新的掩码，只保留满足条件的区域
BW_filtered = false(size(data_filtered));
for i = 1:length(validIdx)
    BW_filtered(CC.PixelIdxList{validIdx(i)}) = true;
end

iso = isosurface(gridX, gridY, gridZ,permute(BW_filtered, [2 1 3]), 0.5);
isosurface_handle = patch(ax,iso, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', alpha);
hold on
% contour(gridX_2D, gridY_2D, dataVolume(:,:,1)', [value_threshold value_threshold])

hold on


% lightangle(45,-45);
% lighting('gouraud');

set(gca, "ZDir", "reverse");
xlabel("longitude");
ylabel("latitude");
zlabel("depth (meter)")
daspect([1,1, 100]);

grid on

end

