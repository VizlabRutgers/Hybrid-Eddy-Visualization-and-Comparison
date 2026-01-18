function [ow_val] = calcOW(xSize, ySize, u_val, v_val)
%CALCOW 此处显示有关此函数的摘要
%   此处显示详细说明
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


    velocity_mag_val = sqrt(u_val.^2 + v_val.^2);
    
    s_n(s_n==0) = NaN;
    s_s(s_n==0) = NaN;
    omega(s_n==0) = NaN;
    ow_val = s_s+s_n-omega;
    ow_val = ow_val./(velocity_mag_val.^2);
end

