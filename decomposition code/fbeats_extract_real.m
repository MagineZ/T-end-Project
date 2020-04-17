function [I2_orig,I2_orig_morph, fbeats_p, fbeats_q, R_amp, S_amp, tfrrF, HR_fe, tfrtic, t] = fbeats_extract(I2_orig,I2_orig_morph,fs,basicTF, advTF, cepR, P, lam_curve, lam_beat, iflog, HR_ma)


I2_orig_len = floor(length(I2_orig)/200)*200;
I2 = resample(I2_orig,100,fs);
basicTF.win = round(basicTF.win*3/5) ;
if iflog
    gg = log(1+abs(I2)) ; gg = gg-mean(gg) ;
else
    gg = I2; gg = gg-mean(gg);
end

%% apply the de-shape on the rough fECG to get the fetal HR
[~, ~, ~, tfrrF, ~, ~, tfrtic, t] = CFPH(gg, basicTF, advTF, cepR, P);



%% supress the possible residure of the mECG
for ti = 1:size(tfrrF, 2)
    idx = [round(HR_ma(ti)*0.94):round(HR_ma(ti)*1.06)] ;
    if max(HR_ma>=size(tfrrF,1))
        continue;
    else
        tfrrF(idx, ti) = tfrrF(idx, ti)/10 ;
    end
    
end

HR_fe = Extract_fhr(tfrrF, basicTF, lam_curve);
HR_fe3 = interp1(200:200:length(I2_orig), HR_fe, 1:length(I2_orig),'pchip','extrap') ;


%% apply the beat tracking to get the fetal HR
flocsf1p = beat_simple(I2, 100, HR_fe3.*basicTF.fr, lam_beat);
flocsf1q = beat_simple(-I2, 100, HR_fe3.*basicTF.fr, lam_beat);


%% get the fetal polarity
if abs(median(I2(flocsf1p))) > abs(median(I2(flocsf1q)))
    Po = 1 ; flocsf = flocsf1p ;
else
    Po = -1 ; flocsf = flocsf1q ;
    fprintf('\t\t*** reverse the fetal pole\n') ;
end

I2 = Po.*I2;
I2_orig = Po.*I2_orig ;
I2_orig_morph = Po.*I2_orig_morph ;

flocsf = RRconstraint(flocsf, I2, 100, 0.25);

I2 = I2_orig;
flocsf = round(flocsf*10/1);
fbeats_p = [];
fbeats_q = [];
SearchLen_p = round(48./mean(HR_fe3*basicTF.fr)) ;
SearchLen_q = 48;

for ii = 1:length(flocsf)
    [~, idx] = max(I2(max([1 flocsf(ii)-SearchLen_p]):min([flocsf(ii)+SearchLen_p length(I2)])));
    fbeats_p(ii) = flocsf(ii)-SearchLen_p+idx-1;
    [~, idx2] = min(I2(max([1 flocsf(ii)-SearchLen_q]):min([flocsf(ii)+SearchLen_q length(I2)])));
    fbeats_q(ii) = flocsf(ii)-SearchLen_q+idx2-1;
end
%fbeats = fbeats(fbeats>0); fbeats = fbeats(fbeats<length(I2_orig));
ind = find(fbeats_p>0 & fbeats_p<I2_orig_len & fbeats_q>0 & fbeats_q<I2_orig_len);
fbeats_p = fbeats_p(ind);
fbeats_q = fbeats_q(ind);
R_amp = median(I2(fbeats_p));
S_amp = median(I2(fbeats_q));

I2_orig = Po.*I2_orig ;
I2_orig_morph = Po.*I2_orig_morph ;


end



