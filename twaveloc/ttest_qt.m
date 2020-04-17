clear all; close all;
load('For_stat.mat');

disp('0 vs 1 ME')
set = 1
[~,p] = ttest(M1Z(:,1),M1Z(:,2),'Tail','right');
disp(['Zhang, set =', int2str(set)',': ' ,num2str(p)])
[~,p] = ttest(M1C(:,1),M1C(:,2),'Tail','right');
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(p)])
[~,p] = ttest(M1W(:,1),M1W(:,2),'Tail','right');
disp(['Martinez, set =' ,int2str(set),': ',num2str(p)])
set = 2
[~,p] = ttest(M2Z(:,1),M2Z(:,2),'Tail','right');
disp(['Zhang, set =', int2str(set)',': ' ,num2str(p)])
[~,p] = ttest(M2C(:,1),M2C(:,2),'Tail','right');
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(p)])
[~,p] = ttest(M2W(:,1),M2W(:,2),'Tail','right');
disp(['Martinez, set =' ,int2str(set),': ',num2str(p)])
set = 3
[~,p] = ttest(M3Z(:,1),M3Z(:,2),'Tail','right');
disp(['Zhang, set =', int2str(set)',': ' ,num2str(p)])
[~,p] = ttest(M3C(:,1),M3C(:,2),'Tail','right');
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(p)])
[~,p] = ttest(M3W(:,1),M3W(:,2),'Tail','right');
disp(['Martinez, set =' ,int2str(set),': ',num2str(p)])
%{
disp('0 vs 1 SD')
set = 1
disp(['Zhang, set =', int2str(set)',': ' ,num2str(ttest(S1Z(:,1),S1Z(:,2)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(ttest(S1C(:,1),S1C(:,2)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(ttest(S1W(:,1),S1W(:,2)))])
set = 2
disp(['Zhang, set =', int2str(set)',': ' ,num2str(ttest(S2Z(:,1),S2Z(:,2)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(ttest(S2C(:,1),S2C(:,2)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(ttest(S2W(:,1),S2W(:,2)))])
set = 3
disp(['Zhang, set =', int2str(set)',': ' ,num2str(ttest(S3Z(:,1),S3Z(:,2)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(ttest(S3C(:,1),S3C(:,2)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(ttest(S3W(:,1),S3W(:,2)))])
%}
disp('0 vs 0.85 ME')
set = 1
[~,p] = ttest(M1Z(:,1),M1Z(:,3),'Tail','right');
disp(['Zhang, set =', int2str(set)',': ' ,num2str(p)])
[~,p] = ttest(M1C(:,1),M1C(:,3),'Tail','right');
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(p)])
[~,p] = ttest(M1W(:,1),M1W(:,3),'Tail','right');
disp(['Martinez, set =' ,int2str(set),': ',num2str(p)])
set = 2
[~,p] = ttest(M2Z(:,1),M2Z(:,3),'Tail','right');
disp(['Zhang, set =', int2str(set)',': ' ,num2str(p)])
[~,p] = ttest(M2C(:,1),M2C(:,3),'Tail','right');
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(p)])
[~,p] = ttest(M2W(:,1),M2W(:,3),'Tail','right');
disp(['Martinez, set =' ,int2str(set),': ',num2str(p)])
set = 3
[~,p] = ttest(M3Z(:,1),M3Z(:,3),'Tail','right');
disp(['Zhang, set =', int2str(set)',': ' ,num2str(p)])
[~,p] = ttest(M3C(:,1),M3C(:,3),'Tail','right');
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(p)])
[~,p] = ttest(M3W(:,1),M3W(:,3),'Tail','right');
disp(['Martinez, set =' ,int2str(set),': ',num2str(p)])
%{
disp('0 vs 0.85 SD')
set = 1
disp(['Zhang, set =', int2str(set)',': ' ,num2str(ttest(S1Z(:,1),S1Z(:,3)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(ttest(S1C(:,1),S1C(:,3)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(ttest(S1W(:,1),S1W(:,3)))])
set = 2
disp(['Zhang, set =', int2str(set)',': ' ,num2str(ttest(S2Z(:,1),S2Z(:,3)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(ttest(S2C(:,1),S2C(:,3)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(ttest(S2W(:,1),S2W(:,3)))])
set = 3
disp(['Zhang, set =', int2str(set)',': ' ,num2str(ttest(S3Z(:,1),S3Z(:,3)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(ttest(S3C(:,1),S3C(:,3)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(ttest(S3W(:,1),S3W(:,3)))])
%}



disp('1 vs 0.85')
set = 1
[~,p] = ttest(M1Z(:,2),M1Z(:,3),'Tail','right');
disp(['Zhang, set =', int2str(set)',': ' ,num2str(p)])
[~,p] = ttest(M1C(:,2),M1C(:,3),'Tail','right');
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(p)])
[~,p] = ttest(M1W(:,2),M1W(:,3),'Tail','right');
disp(['Martinez, set =' ,int2str(set),': ',num2str(p)])
set = 2
[~,p] = ttest(M2Z(:,2),M2Z(:,3),'Tail','right');
disp(['Zhang, set =', int2str(set)',': ' ,num2str(p)])
[~,p] = ttest(M2C(:,2),M2C(:,3),'Tail','right');
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(p)])
[~,p] = ttest(M2W(:,2),M2W(:,3),'Tail','right');
disp(['Martinez, set =' ,int2str(set),': ',num2str(p)])
set = 3
[~,p] = ttest(M3Z(:,2),M3Z(:,3),'Tail','right');
disp(['Zhang, set =', int2str(set)',': ' ,num2str(p)])
[~,p] = ttest(M3C(:,2),M3C(:,3),'Tail','right');
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(p)])
[~,p] = ttest(M3W(:,2),M3W(:,3),'Tail','right');
disp(['Martinez, set =' ,int2str(set),': ',num2str(p)])
%{
disp('1 vs 0.85 SD')
set = 1
disp(['Zhang, set =', int2str(set)',': ' ,num2str(ttest(S1Z(:,2),S1Z(:,3)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(ttest(S1C(:,2),S1C(:,3)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(ttest(S1W(:,2),S1W(:,3)))])
set = 2
disp(['Zhang, set =', int2str(set)',': ' ,num2str(ttest(S2Z(:,2),S2Z(:,3)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(ttest(S2C(:,2),S2C(:,3)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(ttest(S2W(:,2),S2W(:,3)))])
set = 3
disp(['Zhang, set =', int2str(set)',': ' ,num2str(ttest(S3Z(:,2),S3Z(:,3)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(ttest(S3C(:,2),S3C(:,3)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(ttest(S3W(:,2),S3W(:,3)))])

%}