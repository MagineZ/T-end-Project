%clear all;
%close all;
addpath(genpath('D:\duke box\fECG_project\Code Archive3\Code Archive3'))
addpath(genpath('D:\duke box\fECG_project\ecg-kit-0.1.6\ecg-kit-0.1.6\common'));
addpath(genpath('D:\duke box\fECG_project\MFEToolbox'))
noise_mth = "gaussian";
%noise_mth = "arma11";
sqi_th=0.9;
rng(102);
Table_g2 = [];
Table_g3 = [];
for snr = [10,5]
    meZ = []; sdZ = []; meZ_sh = []; sdZ_sh = [];
    meC = []; sdC = []; meC_sh = []; sdC_sh = [];
    meW = []; sdW = []; meW_sh = []; sdW_sh = [];
    for i = 1:20
        [a,b,c,d,e,f,g,h] = blpsample2_all_noise0405(snr,noise_mth,"Zhang",sqi_th);
       meZ = [meZ;[a,e]]; sdZ = [sdZ;[b,f]];  meZ_sh = [meZ_sh;[c,g]]; sdZ_sh = [sdZ_sh;[d,h]];
       [a,b,c,d,e,f,g,h] = blpsample2_all_noise0405(snr,noise_mth,"Carlos",sqi_th);
       meC = [meC;[a,e]]; sdC = [sdC;[b,f]];  meC_sh = [meC_sh;[c,g]]; sdC_sh = [sdC_sh;[d,h]];
       [a,b,c,d,e,f,g,h] = blpsample2_all_noise0405(snr,noise_mth,"Wavedet",sqi_th);
       meW = [meW;[a,e]]; sdW = [sdW;[b,f]];  meW_sh = [meW_sh;[c,g]]; sdW_sh = [sdW_sh;[d,h]];
    end
    save(['loop_test2_gaussian' int2str(snr) 'dB.mat'])
    Table_g2 = [Table_g2;nanmean(meZ(:,1)) ,nanmean(sdZ(:,1)), nanmean(meZ_sh(:,1)) ,nanmean(sdZ_sh(:,1)), nanmean(meC(:,1)), nanmean(sdC(:,1)), nanmean(meC_sh(:,1)), nanmean(sdC_sh(:,1)),nanmean(meW(:,1)), nanmean(sdW(:,1)),nanmean(meW_sh(:,1)), nanmean(sdW_sh(:,1))];
    Table_g3 = [Table_g3;nanmean(meZ(:,2)) ,nanmean(sdZ(:,2)), nanmean(meZ_sh(:,2)) ,nanmean(sdZ_sh(:,2)), nanmean(meC(:,2)), nanmean(sdC(:,2)), nanmean(meC_sh(:,2)), nanmean(sdC_sh(:,2)),nanmean(meW(:,2)), nanmean(sdW(:,2)),nanmean(meW_sh(:,2)), nanmean(sdW_sh(:,2))];
end
