function   fbeats = RR_artifact_detect_fusion(fbeats, median_rr,th3_1,th3_2)
RR = diff(fbeats);
i = 2;
while i <= length(RR)
    
   
    median_rr = round(median_rr);
    if RR(i)-median_rr <th3_1
              
        fbeats(i+1) = [];
    
    elseif RR(i)-median_rr >th3_2
           
        fbeats = [fbeats(1:i),fbeats(i)+median_rr,fbeats(i+1:end)];
      
    else
         i=i+1;       
    end
    
    RR = fbeats(2:end)-fbeats(1:end-1);
   
end