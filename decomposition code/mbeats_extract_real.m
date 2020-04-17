function [x0, x0_real,mbeats_p, mbeats_q, R_amp, S_amp, tfrrM, HR_ma, tfrtic, t] = mbeats_extract_real(x0,x0_real,fs,basicTF, advTF, cepR, P, lam_curve, lam_beat, iflog)
x0_len = floor(length(x0)/200)*200;
x0 = x0(1:x0_len);

x1 = resample(x0,100,fs);
if iflog
    gg = log(1+abs(x1)) ; gg = gg-mean(gg) ;
else
    gg = x1; gg = gg-mean(gg);
end

[~, ~, ~, tfrrM, ~, ~, tfrtic, t] = CFPH(gg, basicTF, advTF, cepR, P);

HR_ma = Extract_mhr(tfrrM, basicTF, lam_curve);
%HR_ma2 = interp1(200:200:60000, HR_ma, 1:60000,'pchip','extrap') ;
HR_ma2 = interp1(200:200:x0_len, HR_ma, 1:x0_len,'pchip','extrap') ;
%% Get maternal R peaks by DP
%mlocsf1 = beat_simple(x1, 100, HR_ma2.*basicTF.fr, lam_beat);

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
x0_real = Po.*x0_real;
%mlocsf1  = beat_simple3(x1, 100, HR_ma2.*basicTF.fr, lam_beat,mlocsf1, 40);%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mlocsf1 = RRconstraint(mlocsf1, x1, 100, 0.25);

%% Get maternal ECG waveforms from nonlocal median
x_up = x0;
mlocsf = round(mlocsf1*10);
mbeats_p = [];
mbeats_q = [];
SearchLen_p = round(48./mean(HR_ma2*basicTF.fr)) ;
SearchLen_q = 48;
for ii = 1:length(mlocsf)
    [~, idx] = max(x_up(max([mlocsf(ii)-SearchLen_p 1]):min([mlocsf(ii)+SearchLen_p length(x_up)]))) ;
    mbeats_p(ii) = mlocsf(ii)-SearchLen_p+idx-1 ;
    [~, idx2] = min(x_up(max([mlocsf(ii)-SearchLen_q 1]):min([mlocsf(ii)+SearchLen_q length(x_up)]))) ;
    mbeats_q(ii) = mlocsf(ii)-SearchLen_q+idx2-1 ;
end
ind = find(mbeats_p>0 & mbeats_p<x0_len & mbeats_q>0 & mbeats_q<x0_len);
mbeats_p = mbeats_p(ind);
mbeats_q = mbeats_q(ind);
R_amp = median(x_up(mbeats_p));
S_amp = median(x_up(mbeats_q));

x0 = Po.*x0 ;
x0_real = Po.*x0_real;

