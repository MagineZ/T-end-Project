clear all; close all;
noise = 'gaussian'
snr = 5
load(['loop_test_' noise int2str(snr) 'dB.mat'])
Table1 = [];
Table1 = [Table1, median(meZ), mad(meZ), median(meZ_sh), mad(meZ_sh)];
Table1 = [Table1, median(meC), mad(meC), median(meC_sh), mad(meC_sh)];
Table1 = [Table1, median(meW), mad(meW), median(meW_sh), mad(meW_sh)];
Table1q = [];
p1 = 0.025; p2 = 0.975;
Table1q = [Table1q, quantile(meZ,p1), quantile(meZ,p2), quantile(meZ_sh,p1), quantile(meZ_sh,p2)];
Table1q = [Table1q, quantile(meC,p1), quantile(meC,p2), quantile(meC_sh,p1), quantile(meC_sh,p2)];
Table1q = [Table1q, quantile(meW,p1), quantile(meW,p2), quantile(meW_sh,p1), quantile(meW_sh,p2)];

disp('0 vs 0.85')
disp(['Zhang, set =1: ' ,num2str(signrank(meZ,meZ_sh))])
disp(['Carlos, set =1: ' ,num2str(signrank(meC,meC_sh))])
disp(['Martinez, set =1: ',num2str(signrank(meW,meW_sh))])


load(['loop_test2_' noise int2str(snr) 'dB.mat'])
Table2 = [];
Table2 = [Table2, median(meZ(:,1)), mad(meZ(:,1)), median(meZ_sh(:,1)), mad(meZ_sh(:,1))];
Table2 = [Table2, median(meC(:,1)), mad(meC(:,1)), median(meC_sh(:,1)), mad(meC_sh(:,1))];
Table2 = [Table2, median(meW(:,1)), mad(meW(:,1)), median(meW_sh(:,1)), mad(meW_sh(:,1))];
Table3 = [];
Table3 = [Table3, median(meZ(:,2)), mad(meZ(:,2)), median(meZ_sh(:,2)), mad(meZ_sh(:,2))];
Table3 = [Table3, median(meC(:,2)), mad(meC(:,2)), median(meC_sh(:,2)), mad(meC_sh(:,2))];
Table3 = [Table3, median(meW(:,2)), mad(meW(:,2)), median(meW_sh(:,2)), mad(meW_sh(:,2))];
Table2q = [];
Table2q = [Table2q, quantile(meZ(:,1),p1), quantile(meZ(:,1),p2), quantile(meZ_sh(:,1),p1), quantile(meZ_sh(:,1),p2)];
Table2q = [Table2q, quantile(meC(:,1),p1), quantile(meC(:,1),p2), quantile(meC_sh(:,1),p1), quantile(meC_sh(:,1),p2)];
Table2q = [Table2q, quantile(meW(:,1),p1), quantile(meW(:,1),p2), quantile(meW_sh(:,1),p1), quantile(meW_sh(:,1),p2)];
Table3q = [];
Table3q = [Table3q, quantile(meZ(:,2),p1), quantile(meZ(:,2),p2), quantile(meZ_sh(:,2),p1), quantile(meZ_sh(:,2),p2)];
Table3q = [Table3q, quantile(meC(:,2),p1), quantile(meC(:,2),p2), quantile(meC_sh(:,2),p1), quantile(meC_sh(:,2),p2)];
Table3q = [Table3q, quantile(meW(:,2),p1), quantile(meW(:,2),p2), quantile(meW_sh(:,2),p1), quantile(meW_sh(:,2),p2)];


load(['loop_test2_' noise int2str(snr) 'dB.mat'])
disp(['Zhang, set =2: ' ,num2str(signrank(meZ(:,1),meZ_sh(:,1)))])
disp(['Carlos, set =2: ' ,num2str(signrank(meC(:,1),meC_sh(:,1)))])
disp(['Martinez, set =2: ',num2str(signrank(meW(:,1),meW_sh(:,1)))])

disp(['Zhang, set =3: ' ,num2str(signrank(meZ(:,2),meZ_sh(:,2)))])
disp(['Carlos, set =3: ' ,num2str(signrank(meC(:,2),meC_sh(:,2)))])
disp(['Martinez, set =3: ',num2str(signrank(meW(:,2),meW_sh(:,2)))])
