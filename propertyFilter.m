function [outputEddy] = propertyFilter(inputEddy,spaceLimit, radiusLimit, streamRidge, x_val, y_val, northFlag)
%PROPERTYFILTER 此处显示有关此函数的摘要
%   此处显示详细说明
inputEddy(inputEddy(:,1)<spaceLimit.x0,:) = [];
inputEddy(inputEddy(:,1)>spaceLimit.x1,:) = [];
inputEddy(inputEddy(:,2)<spaceLimit.y0,:) = [];
inputEddy(inputEddy(:,2)>spaceLimit.y1,:) = [];

inputEddy(inputEddy(:,13)<radiusLimit.lower,:) = [];
inputEddy(inputEddy(:,13)>radiusLimit.upper,:) = [];


[~, lower_x] = min(abs(x_val - spaceLimit.x0));
[~, upper_x] = min(abs(x_val - spaceLimit.x1));
[~, lower_y] = min(abs(y_val - spaceLimit.y0));
[~, upper_y] = min(abs(y_val - spaceLimit.y1));


eddyX = inputEddy(:,1);
eddyY = inputEddy(:,2);
eddyFrame = inputEddy(:,15);


[~, loc_x] = ismember(round(eddyX,4), round(x_val,4));
[~, loc_y] = ismember(round(eddyY,4), round(y_val,4));

loc_x_correlated = loc_x-lower_x;
loc_y_correlated = loc_y-lower_y;


outputEddy_filterFlag = ones(size(inputEddy,1),1);

if(streamRidge~=0)
    for index = 1:1:size(inputEddy,1)
        if(northFlag == 1)
            if(streamRidge(loc_x_correlated(index), eddyFrame(index))<loc_y_correlated(index))
                outputEddy_filterFlag(index)=0;
            end
        elseif(northFlag == 0)
            if(streamRidge(loc_x_correlated(index), eddyFrame(index))>loc_y_correlated(index))
                outputEddy_filterFlag(index)=0;
            end
        end
    end
end



outputEddy = inputEddy(outputEddy_filterFlag==1,:);



end

