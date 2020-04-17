function x0 = ECG_detrend(x0, num_med_s, num_med_l, if_morph, if_smooth)
%detrend
%if_morph = 1 : preserve morphology
%         = 0 : otherwise
if nargin < 5
    if_smooth = 1;
end
    
if if_morph
    x0_trend = movmedian(x0,num_med_s);% medfilt2(x0, [1, num_med_s]);
    x0_trend = movmedian(x0_trend,num_med_l);% medfilt2(x0_trend, [1, num_med_l]);
   
   
else
    x0_trend = movmedian(x0,num_med_s);% medfilt2(x0, [1, num_med_s]);
end

if if_smooth
    x0_trend = smooth(x0_trend,10,'loess')';
end
 x0 = x0 - x0_trend;
end