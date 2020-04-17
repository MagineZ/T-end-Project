function Tend = detect_Tend_Carlos(ECG,Tpeaks,fs)
%Twin = round(fs * 0.2);
Twin = round(median(diff(Tpeaks)) * 0.2);
%Twin = [Twin(1),Twin];
x0 = smooth(ECG,fs*0.03);
%x0 = ECG;
Tend = [];
for i = 1:length(Tpeaks)

    if Tpeaks(i)+Twin-1 > length(x0)
        tw = x0(Tpeaks(i):end);
    else
        tw = x0(Tpeaks(i):(Tpeaks(i)+Twin-1));
    end
    tw_diff = abs(diff(tw));
    [~,t_m] = max(tw_diff);
    y_m = tw(t_m);
    t_r = length(tw); 
    A = 0;
    tend_loc = 1;
    for j = t_m:t_r 
        current_A =  abs(0.5*(y_m - tw(j))*(2*t_r - j - t_m));
        if current_A >A
            A = current_A;
            tend_loc = j;
        end
    end
    Tend = [Tend,Tpeaks(i)+tend_loc-1];
    
end
end