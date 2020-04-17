function   fbeats = RR_artifact_detect(fbeats, num_nb)
RR = diff(fbeats);
i = 2;
while i <= length(RR)
    
    if i <= ceil(num_nb/2)
        median_rr = median(RR(1:num_nb));
    elseif i>= length(RR) - floor(num_nb/2)
        median_rr = median(RR(end-num_nb+1:end));
    else
        median_rr = median(RR(i-floor(num_nb/2):i+floor(num_nb/2)));
    end
    median_rr = round(median_rr);
    if (abs(RR(i)-median_rr)/median_rr>0.2) && (abs(RR(i)-RR(i-1))/RR(i-1)>0.2)
        if RR(i) < median_rr
           fbeats(i+1) = [];
        else
           fbeats = [fbeats(1:i),fbeats(i)+median_rr,fbeats(i+1:end)];
        end
    else
         i=i+1;       
    end
    RR = fbeats(2:end)-fbeats(1:end-1);
   
   
end