clear all; close all;
noise = 'gaussian'
snr = 10
load(['loop_test_' noise int2str(snr) 'dB.mat'])

disp('0 vs 0.85')
[~,p] = ttest(meZ,meZ_sh,'Tail','right');
disp(['Zhang, set =1: ' ,num2str(p)])
[~,p] = ttest(meC,meC_sh,'Tail','right');
disp(['Carlos, set =1: ' ,num2str(p)])
[~,p] = ttest(meW,meW_sh,'Tail','right');
disp(['Martinez, set =1: ',num2str(p)])

load(['loop_test2_' noise int2str(snr) 'dB.mat'])
[~,p] = ttest(meZ(:,1),meZ_sh(:,1),'Tail','right');
disp(['Zhang, set =2: ' ,num2str(p)])
[~,p] = ttest(meC(:,1),meC_sh(:,1),'Tail','right');
disp(['Carlos, set =2: ' ,num2str(p)])
[~,p] = ttest(meW(:,1),meW_sh(:,1),'Tail','right');
disp(['Martinez, set =2: ',num2str(p)])

[~,p] = ttest(meZ(:,2),meZ_sh(:,2),'Tail','right');
disp(['Zhang, set =3: ' ,num2str(p)])
[~,p] = ttest(meC(:,2),meC_sh(:,2),'Tail','right');
disp(['Carlos, set =3: ' ,num2str(p)])
[~,p] = ttest(meW(:,2),meW_sh(:,2),'Tail','right');
disp(['Martinez, set =3: ',num2str(p)])
