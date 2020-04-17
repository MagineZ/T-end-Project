function [alpha_mbeats_c,alpha_sig_c,HR_ma] = mbeats_modify_tfrrMmedian(alpha_tfrrM_c, alpha_sig_c,basicTF, lam_curve,lam_beat,fs,iflog)


tfrrM = median(alpha_tfrrM_c,3);
x0_len = floor(length(alpha_sig_c{1})/200)*200;

HR_ma = Extract_mhr(tfrrM, basicTF, lam_curve);
HR_ma2 = interp1(200:200:x0_len, HR_ma, 1:x0_len,'pchip','extrap') ;
[~,n] = size(alpha_sig_c);
alpha_mbeats_c = {};
for i = 1:n
    x0 = alpha_sig_c{i};
    x0 = x0(1:x0_len);
    x1 = resample(x0,100,fs);
    if iflog
        gg = log(1+abs(x1)) ; gg = gg-mean(gg) ;
    else
        gg = x1; gg = gg-mean(gg);
    end
    
    mlocsf1p = beat_simple(x1, 100, HR_ma2.*basicTF.fr, lam_beat);
    mlocsf1q = beat_simple(-x1, 100, HR_ma2.*basicTF.fr, lam_beat);
    %mlocsf1p = RRconstraint(mlocsf1p, x1, 100, 0.25);
    %mlocsf1q = RRconstraint(mlocsf1q, -x1, 100, 0.25);
    
    
    if abs(median(x1(mlocsf1p))) > abs(median(x1(mlocsf1q)))
        %if std(mlocsf1p) < std(mlocsf1q)
        Po = 1 ; mlocsf1 = mlocsf1p ;
    else
        Po = -1 ; mlocsf1 = mlocsf1q ;
    end
    x1 = Po.*x1 ;
    x0 = Po.*x0 ;
    
    mlocsf1 = RRconstraint(mlocsf1, x1, 100, 0.25);
    
    %% Get maternal ECG waveforms from nonlocal median
    x_up = x0;
    mlocsf = round(mlocsf1*10);
    mbeats = [];
    SearchLen_p = round(48./mean(HR_ma2*basicTF.fr)) ;
    for ii = 1:length(mlocsf)
        [~, idx] = max(x_up(max([mlocsf(ii)-SearchLen_p 1]):min([mlocsf(ii)+SearchLen_p length(x_up)]))) ;
        mbeats(ii) = mlocsf(ii)-SearchLen_p+idx-1 ;
    end
    ind = find(mbeats>0 & mbeats<x0_len);
    mbeats = mbeats(ind);
    alpha_mbeats_c{i} = mbeats;
    alpha_sig_c{i} = x0;
end