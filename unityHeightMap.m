function [outputArg1,outputArg2] = unityHeightMap(srcData)
%UNITYHEIGHTMAP 此处显示有关此函数的摘要
%   此处显示详细说明
x_val = double(ncread(srcData, "XC"));
y_val = double(ncread(srcData, "YC"));
z_val = double(ncread(srcData, "Z_MIT40"));
resolution = 0.04;


FrameNum=1;
startLoc = [1,1,1,FrameNum];
count = [length(x_val),length(y_val),length(z_val),1];
stride = [1,1,1,1];



coastRegion = ncread(srcData,"TEMP",startLoc, count, stride);
coastHeightMap = zeros(length(x_val),length(y_val));

[~,coastHeightMap] = min(coastRegion,[],3);
coastHeightMap = imcomplement(im2uint8(mat2gray(coastHeightMap)));
imshow(coastHeightMap);


fid=fopen("coastHeight_RedSea.raw",'w+');
cnt=fwrite(fid,coastHeightMap,'uint8');
fclose(fid);
end

