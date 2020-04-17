%clear all;
%close all;
addpath(genpath('D:\duke box\fECG_project\Code Archive3\Code Archive3'))
addpath(genpath('D:\duke box\fECG_project\ecg-kit-0.1.6\ecg-kit-0.1.6\common'));
addpath(genpath('D:\duke box\fECG_project\MFEToolbox'))
noise_mth = "gaussian";
%noise_mth = "arma11";
sqi_th = 0.9;
rng(100);
Table_g1 = [];
for snr = [10,5]
    meZ = []; sdZ = []; meZ_sh = []; sdZ_sh = [];
    meC = []; sdC = []; meC_sh = []; sdC_sh = [];
    meW = []; sdW = []; meW_sh = []; sdW_sh = [];
    for i = 1:20
       [a,b,c,d]=blpsample1_all_noise0405(snr,noise_mth,"Zhang",sqi_th);
       meZ = [meZ,a]; sdZ = [sdZ,b];  meZ_sh = [meZ_sh,c]; sdZ_sh = [sdZ_sh,d];
       [a,b,c,d]=blpsample1_all_noise0405(snr,noise_mth,"Carlos",sqi_th);
       meC = [meC,a]; sdC = [sdC,b];  meC_sh = [meC_sh,c]; sdC_sh = [sdC_sh,d];
       [a,b,c,d]=blpsample1_all_noise0405(snr,noise_mth,"Wavedet",sqi_th);
       meW = [meW,a]; sdW = [sdW,b];  meW_sh = [meW_sh,c]; sdW_sh = [sdW_sh,d];
        
    end
    save(['loop_test_gaussian' int2str(snr) 'dB.mat'])
    Table_g1 = [Table_g1;nanmean(meZ) ,nanmean(sdZ), nanmean(meZ_sh) ,nanmean(sdZ_sh), nanmean(meC), nanmean(sdC), nanmean(meC_sh), nanmean(sdC_sh),nanmean(meW), nanmean(sdW),nanmean(meW_sh), nanmean(sdW_sh)];
end
