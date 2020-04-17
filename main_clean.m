clear;
addpath(genpath(pwd))
q = 0.9;
fs =1000;
[uerr1,uerr2,bsqi_Tend_pb_set1,tsqi_set1,bsqi_set1]  = blpsample1_all_sqi_0406('Zhang',q,fs);
M1Z = [uerr1,uerr2];

[uerr1,uerr2]  = blpsample1_all_sqi_0406('Carlos',q,fs);
M1C = [uerr1,uerr2];

[uerr1,uerr2]  = blpsample1_all_sqi_0406('Wavedet',q,fs);
M1W = [uerr1,uerr2];

[uerr1,uerr1_sh,uerr2,uerr2_sh,bsqi_Tend_pb_set2,bsqi_Tend_pb_set3,tsqi_set2,tsqi_set3,bsqi_set2,bsqi_set3] = blpsample2_all_sqi_0406('Zhang',q,fs);
M2Z = [uerr1,uerr1_sh];
M3Z = [uerr2,uerr2_sh];


[uerr1,uerr1_sh,uerr2,uerr2_sh] = blpsample2_all_sqi_0406('Carlos',q,fs);

M2C = [uerr1,uerr1_sh];
M3C = [uerr2,uerr2_sh];


[uerr1,uerr1_sh,uerr2,uerr2_sh] = blpsample2_all_sqi_0406('Wavedet',q,fs);
M2W = [uerr1,uerr1_sh];
M3W = [uerr2,uerr2_sh];ts

save(['QT_stat_',int2str(fs)],'M1Z','M1C','M1W','M2Z','M2C','M2W','M3Z','M3C','M3W','bsqi_Tend_pb_set1','bsqi_Tend_pb_set2','bsqi_Tend_pb_set3','tsqi_set1','tsqi_set2','tsqi_set3','bsqi_set1','bsqi_set2','bsqi_set3')