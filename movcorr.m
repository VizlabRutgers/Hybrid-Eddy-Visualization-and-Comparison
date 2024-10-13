function [avgCorr] = movcorr(x,y,n)
    
    if(n>length(x))
        error("The sliding window is supposed to be larger than the array");
    end

    for i = 1:1:length(x)
    
        lower_bound = i-(n-1)/2;
        upper_bound = i+(n-1)/2;
    
        if(lower_bound>0 && upper_bound<=length(x))
            sample_x = x(lower_bound : upper_bound);
            sample_y = y(lower_bound : upper_bound);
        elseif(lower_bound<1)
            
    
            sample_x = [zeros(1,1-lower_bound),x(1:upper_bound)];
            sample_y = [zeros(1,1-lower_bound),y(1:upper_bound)];
        elseif(upper_bound>length(x))
            sample_x = [x(lower_bound:length(x)),zeros(1,upper_bound-length(x))];
            sample_y = [y(lower_bound:length(x)),zeros(1,upper_bound-length(x))];
        end
    
    
    
        corr(i) = fixedCorr(sample_x,sample_y);
    

    end
    avgCorr = mean(corr((1+(n-1)/2):(length(corr)-(n-1)/2)));
end



