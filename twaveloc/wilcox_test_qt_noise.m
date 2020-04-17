clear all; close all;
noise = 'arma11'
snr = 5
load(['loop_test_' noise int2str(snr) 'dB.mat'])

disp('0 vs 0.85')
disp(['Zhang, set =1: ' ,num2str(signrank(meZ,meZ_sh))])
disp(['Carlos, set =1: ' ,num2str(signrank(meC,meC_sh))])
disp(['Martinez, set =1: ',num2str(signrank(meW,meW_sh))])

load(['loop_test2_' noise int2str(snr) 'dB.mat'])
disp(['Zhang, set =2: ' ,num2str(signrank(meZ(:,1),meZ_sh(:,1)))])
disp(['Carlos, set =2: ' ,num2str(signrank(meC(:,1),meC_sh(:,1)))])
disp(['Martinez, set =2: ',num2str(signrank(meW(:,1),meW_sh(:,1)))])

disp(['Zhang, set =3: ' ,num2str(signrank(meZ(:,2),meZ_sh(:,2)))])
disp(['Carlos, set =3: ' ,num2str(signrank(meC(:,2),meC_sh(:,2)))])
disp(['Martinez, set =3: ',num2str(signrank(meW(:,2),meW_sh(:,2)))])
