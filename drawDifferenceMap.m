function [outputArg1,outputArg2] = drawDifferenceMap(statMatrix, srcData, property)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
x_val = double(ncread(srcData, property.x));
y_val = double(ncread(srcData, property.y));
z_val = double(ncread(srcData, property.z));
totalFrames = double(ncread(srcData, property.time));
[x_grid, y_grid] = meshgrid(x_val,y_val);


%% background
f1 = figure();
f1.WindowState = 'maximized';
ax1 = axes(f1);

frame = min(statMatrix(:,5));

startLoc = [1,1,1,frame];
count = [length(x_val),length(y_val),1,1];
stride = [1,1,1,1];
u_val = ncread(srcData, property.u, startLoc, count, stride);
v_val = ncread(srcData, property.v,startLoc, count, stride);

% w_val = ncread(srcpath, "W",startLoc, count, stride);
% eta_val = ncread(srcData, property.eta,[1,1,1], [length(x_val),length(y_val),length(totalFrames)], [1,1,1]);
% temp_val = ncread(srcpath, "TEMP",startLoc, count, stride);
% salinity_val = ncread(srcpath, "SALT",startLoc, count, stride);
vel_mag = sqrt(u_val.^2+v_val.^2);

xSize = length(x_val);
ySize = length(y_val);

u_val_test = u_val(1:xSize-1,:,:) - u_val(2:xSize,:,:);
u_x = u_val;
u_x(1:xSize-1,:,:) = u_val_test;
u_val_test = u_val(:,1:ySize-1,:) - u_val(:,2:ySize,:);
u_y = u_val;
u_y(:,1:ySize-1,:) = u_val_test;

v_val_test = v_val(1:xSize-1,:,:) - v_val(2:xSize,:,:);
v_x = v_val;
v_x(1:xSize-1,:,:) = v_val_test;
v_val_test = v_val(:,1:ySize-1,:) - v_val(:,2:ySize,:);
v_y = v_val;
v_y(:,1:ySize-1,:) = v_val_test;

s_n = (u_x-v_y).^2;
s_s = (u_y+v_x).^2;
omega = (v_x-u_y).^2;


s_n(s_n==0) = NaN;
s_s(s_n==0) = NaN;
omega(s_n==0) = NaN;
ow_val = s_s+s_n-omega;
% Normalization
ow_val = ow_val./(vel_mag.^2);



x_test = x_val(1:length(x_val));
y_test = y_val(1:length(y_val));

u_test = u_val(1:length(x_val),1:length(y_val),1);
v_test = v_val(1:length(x_val),1:length(y_val),1);

u_test(isnan(u_test)) = 0;
v_test(isnan(v_test)) = 0;

vel_mag_test = sqrt(u_test.^2+v_test.^2);



[x_Rgrid2D,y_Rgrid2D] = ndgrid(x_test,y_test);
DataSize_2D = [length(x_val), length(y_val)];

cla(ax1);

%     I0 = uimagesc(ax1, x_test, y_test, (~(vel_mag_test')));
background = pcolor(ax1, x_grid, y_grid, double(vel_mag_test'==0));
background.FaceColor=[0.3 0.3 0.3];
background.EdgeColor = 'none';
alpha(background, double(vel_mag_test'==0)); 

hold on


% color quiver
q = quiver(ax1,x_Rgrid2D(1:3:end,1:3:end),y_Rgrid2D(1:3:end,1:3:end),u_test(1:3:end,1:3:end),v_test(1:3:end,1:3:end), "r", "AutoScaleFactor",1.5,"LineWidth",1);
set(ax1,'YDir','normal');
hold on

x1 = xlabel(ax1,'Logitude');
y1 = ylabel(ax1,'Latitude');
daspect(ax1,[1,1,1]);

ax1.FontSize=24;
x1.FontSize = 48;
y1.FontSize = 48;


%% data


% statMatrixFiltered = statMatrix(statMatrix(:,5) <=10,:);
statMatrixFiltered = statMatrix;

dataSize = statMatrixFiltered(:,1)./max(statMatrixFiltered(:,1));
dataSize(dataSize<0.1) = 0.1;
dataSize = round(dataSize*200);

hold on
scatter(statMatrixFiltered(:,6), statMatrixFiltered(:,7), dataSize, statMatrixFiltered(:,2),'filled');
cb1 = colorbar();
caxis(ax1,[0,0.6]);
cb1.Label.String = "Center Distance ratio";
cb1.TickLabels = strcat(string(cb1.Ticks * 100), '%');

hold off

f2 = figure();
f2.WindowState = 'maximized';
ax2 = axes(f2);



background = pcolor(ax2, x_grid, y_grid, double(vel_mag_test'==0));
background.FaceColor=[0.3 0.3 0.3];
background.EdgeColor = 'none';
alpha(background, double(vel_mag_test'==0)); 

hold on


% color quiver
q = quiver(ax2,x_Rgrid2D(1:3:end,1:3:end),y_Rgrid2D(1:3:end,1:3:end),u_test(1:3:end,1:3:end),v_test(1:3:end,1:3:end), "r", "AutoScaleFactor",1.5,"LineWidth",1);
set(ax2,'YDir','normal');
hold on

dataSize = statMatrixFiltered(:,3)./max(statMatrixFiltered(:,3));
dataSize(dataSize<0.1) = 0.1;
dataSize = round(dataSize*200);



hold on
scatter(ax2,statMatrixFiltered(:,6), statMatrixFiltered(:,7), dataSize, statMatrixFiltered(:,4),'filled');
cb2 = colorbar();
caxis(ax2,[0,3]);
cb2.Label.String = "Area difference ratio";
cb2.TickLabels = strcat(string(cb2.Ticks * 100), '%');

fullmap = parula(256);    % 默认 parula
idx = round(linspace(1,256,3));
cmap3 = fullmap(idx, :);  % 取蓝—绿—黄三色

colormap(cmap3);


x2 = xlabel(ax2,'Logitude');
y2 = ylabel(ax2,'Latitude');
daspect(ax2,[1,1,1]);

ax2.FontSize=24;
x2.FontSize = 48;
y2.FontSize = 48;


f3 = figure();
ax3 = axes(f3);
plot(ax3, statMatrixFiltered(:,1).*(1+statMatrixFiltered(:,2)), statMatrixFiltered(:,3).*(1+statMatrixFiltered(:,4)), "*");


x3 = xlabel(ax3,'Center Distance');
y3 = ylabel(ax3,'Area Difference');

ax3.FontSize=24;
x3.FontSize = 48;
y3.FontSize = 48;



x_origin_idx = interp1(x_val, 1:numel(x_val), statMatrixFiltered(:,6), 'nearest', 'extrap');
y_origin_idx = interp1(y_val, 1:numel(y_val), statMatrixFiltered(:,7), 'nearest', 'extrap');
x_new_idx = interp1(x_val, 1:numel(x_val), statMatrixFiltered(:,8), 'nearest', 'extrap');
y_new_idx = interp1(y_val, 1:numel(y_val), statMatrixFiltered(:,9), 'nearest', 'extrap');


% ssh_origin = eta_val(sub2ind(size(eta_val),x_origin_idx,y_origin_idx,statMatrixFiltered(:,5)));
% ssh_new = eta_val(sub2ind(size(eta_val),x_new_idx,y_new_idx,statMatrixFiltered(:,5)));

% f4 = figure();
% ax4 = axes(f4);
% plot(ax4, ssh_new, statMatrixFiltered(:,2), "*");
% 
% 
% x4 = xlabel(ax4,'ssh');
% y4 = ylabel(ax4,'center Distance Ratio');
% 
% ax4.FontSize=24;
% ylim([0 5])
% x4.FontSize = 48;
% y4.FontSize = 48;


f5 = figure();
ax5 = axes(f5);
% plot(ax5, statMatrixFiltered(:,1).*(1+statMatrixFiltered(:,2)), statMatrixFiltered(:,3).*(1+statMatrixFiltered(:,4)), "*");
plot(ax5, statMatrixFiltered(:,2), statMatrixFiltered(:,4), "*");


x5 = xlabel(ax5,'Center Offset Ratio');
y5 = ylabel(ax5,'Area Difference Ratio');

ax5.FontSize=24;
x5.FontSize = 48;
y5.FontSize = 48;
ylim([0,5]);
end

