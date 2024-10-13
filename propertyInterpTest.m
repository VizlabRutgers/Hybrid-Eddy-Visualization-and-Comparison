clc;
clear all;

% load wind
% figure,
% [sx,sy,sz] = meshgrid(80,20:10:50,0:5:15);
% % quiver3(x,y,z,u,v,w);
% streamribbon(x,y,z,u,v,w,sx,sy,sz);
% view(3);
% shading interp;
% camlight;
% lighting gouraud;

%-------Reading Data---------

Path = "../FT_result/1/Seperated Structures/";
file = dir(fullfile(Path,'*.uocd'));
filenames = {file.name};



startLoc = [1,1,1,1];
count = [500,500,50,1];
stride = [1,1,1,1];


x_val = double(ncread("/home/weiping/data/SciViz/SciVisContest2020/ensembles/0001/COMBINED_2011013100.nc", "XC",1, 500, 1));
y_val = double(ncread("/home/weiping/data/SciViz/SciVisContest2020/ensembles/0001/COMBINED_2011013100.nc", "YC",1, 500, 1));
z_val = double(ncread("/home/weiping/data/SciViz/SciVisContest2020/ensembles/0001/COMBINED_2011013100.nc", "Z_MIT40",1, 50, 1));
u_val = ncread("/home/weiping/data/SciViz/SciVisContest2020/ensembles/0001/COMBINED_2011013100.nc", "U", startLoc, count, stride);
v_val = ncread("/home/weiping/data/SciViz/SciVisContest2020/ensembles/0001/COMBINED_2011013100.nc", "V",startLoc, count, stride);
w_val = ncread("/home/weiping/data/SciViz/SciVisContest2020/ensembles/0001/COMBINED_2011013100.nc", "W",startLoc, count, stride);
eta_val = ncread("/home/weiping/data/SciViz/SciVisContest2020/ensembles/0001/COMBINED_2011013100.nc", "ETA",[1,1,1], [500,500,1], [1,1,1]);
temp_val = ncread("/home/weiping/data/SciViz/SciVisContest2020/ensembles/0001/COMBINED_2011013100.nc", "TEMP",startLoc, count, stride);
salinity_val = ncread("/home/weiping/data/SciViz/SciVisContest2020/ensembles/0001/COMBINED_2011013100.nc", "SALT",startLoc, count, stride);

temp_val_1 = ncread("/home/weiping/data/SciViz/SciVisContest2020/ensembles/0001/COMBINED_2011013100.nc", "TEMP",[1,1,1,1], [500,500,1,1],[1,1,1,1]);
temp_val_2 = ncread("/home/weiping/data/SciViz/SciVisContest2020/ensembles/0001/COMBINED_2011013100.nc", "TEMP",[1,1,1,2], [500,500,1,1],[1,1,1,1]);
temp_val_3 = ncread("/home/weiping/data/SciViz/SciVisContest2020/ensembles/0001/COMBINED_2011013100.nc", "TEMP",[1,1,1,3], [500,500,1,1],[1,1,1,1]);

salt_val_1 = ncread("/home/weiping/data/SciViz/SciVisContest2020/ensembles/0001/COMBINED_2011013100.nc", "SALT",[1,1,1,1], [500,500,1,1],[1,1,1,1]);
salt_val_2 = ncread("/home/weiping/data/SciViz/SciVisContest2020/ensembles/0001/COMBINED_2011013100.nc", "SALT",[1,1,1,2], [500,500,1,1],[1,1,1,1]);
salt_val_3 = ncread("/home/weiping/data/SciViz/SciVisContest2020/ensembles/0001/COMBINED_2011013100.nc", "SALT",[1,1,1,3], [500,500,1,1],[1,1,1,1]);


temp_val_1(temp_val_1==0) = nan;
temp_val_2(temp_val_2==0) = nan;
temp_val_3(temp_val_3==0) = nan;


figure,imagesc(temp_val_1);
figure,imagesc(temp_val_2);
figure,imagesc(temp_val_3);


figure,quiver(u_val(:,:,1)', v_val(:,:,1)');


temp_interp_spatial = temp_val_1/2+temp_val_3/2;

figure,imagesc(temp_interp_spatial);
colorbar();

figure,imagesc(temp_interp_spatial - temp_val_2);
colorbar();
title("difference of averge interpolation and ground truth");


temp_1_fft = fft2(temp_val_1);
temp_3_fft = fft2(temp_val_3);

temp_interp_fft = temp_1_fft/3+2*temp_3_fft/3;
temp_interp_ifft = ifft2(temp_interp_fft);


figure,imagesc(temp_interp_ifft - temp_val_2);
colorbar();


temp_1_surface = temp_val_1(:,:,1);
temp_2_surface = temp_val_2(:,:,1);
temp_3_surface = temp_val_3(:,:,1);


cvx_begin
    variable lambda_temp
    minimize( norm(temp_1_surface*lambda_temp+temp_3_surface*(1-lambda_temp)-temp_2_surface) )
cvx_end

temp_interp_optimization = temp_val_1*lambda_temp+temp_val_3*(1-lambda_temp);

figure,imagesc(temp_interp_optimization);

figure,imagesc(temp_interp_optimization - temp_val_2);
colorbar();
caxis([-1, 1]);
title("difference of optimization interpolation and ground truth");


salt_1_surface = salt_val_1(:,:,1);
salt_2_surface = salt_val_2(:,:,1);
salt_3_surface = salt_val_3(:,:,1);

cvx_begin
    variable lambda_salt
    minimize( norm(salt_1_surface*lambda_salt+salt_3_surface*(1-lambda_salt)-salt_2_surface) )
cvx_end







figure,imagesc(temp_val_1 - temp_1_ifft);
