clear all; close all;
load('For_stat.mat');

M1Z = M1Z*4;M2Z = M2Z*4;M3Z = M3Z*4;
M1C = M1C*4;M2C = M2C*4;M3C = M3C*4;
M1W = M1W*4;M2W = M2W*4;M3W = M3W*4;
Table1 = [];
Table1 = [Table1, median(M1Z(:,1)), mad(M1Z(:,1)),median(M1Z(:,2)), mad(M1Z(:,2)),median(M1Z(:,3)), mad(M1Z(:,3))];
Table1 = [Table1, median(M1C(:,1)), mad(M1C(:,1)),median(M1C(:,2)), mad(M1C(:,2)),median(M1C(:,3)), mad(M1C(:,3))];
Table1 = [Table1, median(M1W(:,1)), mad(M1W(:,1)),median(M1W(:,2)), mad(M1W(:,2)),median(M1W(:,3)), mad(M1W(:,3))];
Table2 = [];
Table2 = [Table2, median(M2Z(:,1)), mad(M2Z(:,1)),median(M2Z(:,2)), mad(M2Z(:,2)),median(M2Z(:,3)), mad(M2Z(:,3))];
Table2 = [Table2, median(M2C(:,1)), mad(M2C(:,1)),median(M2C(:,2)), mad(M2C(:,2)),median(M2C(:,3)), mad(M2C(:,3))];
Table2 = [Table2, median(M2W(:,1)), mad(M2W(:,1)),median(M2W(:,2)), mad(M2W(:,2)),median(M2W(:,3)), mad(M2W(:,3))];
Table3 = [];
Table3 = [Table3, median(M3Z(:,1)), mad(M3Z(:,1)),median(M3Z(:,2)), mad(M3Z(:,2)),median(M3Z(:,3)), mad(M3Z(:,3))];
Table3 = [Table3, median(M3C(:,1)), mad(M3C(:,1)),median(M3C(:,2)), mad(M3C(:,2)),median(M3C(:,3)), mad(M3C(:,3))];
Table3 = [Table3, median(M3W(:,1)), mad(M3W(:,1)),median(M3W(:,2)), mad(M3W(:,2)),median(M3W(:,3)), mad(M3W(:,3))];

Table1q = [];
p1 = 0.025; p2 = 0.975;
Table1q = [Table1q, quantile(M1Z(:,1),p1), quantile(M1Z(:,1),p2), quantile(M1Z(:,2),p1), quantile(M1Z(:,2),p2), quantile(M1Z(:,3),p1), quantile(M1Z(:,3),p2),];
Table1q = [Table1q, quantile(M1C(:,1),p1), quantile(M1C(:,1),p2), quantile(M1C(:,2),p1), quantile(M1C(:,2),p2), quantile(M1C(:,3),p1), quantile(M1C(:,3),p2),];
Table1q = [Table1q, quantile(M1W(:,1),p1), quantile(M1W(:,1),p2), quantile(M1W(:,2),p1), quantile(M1W(:,2),p2), quantile(M1W(:,3),p1), quantile(M1W(:,3),p2),];
Table2q = [];
Table2q = [Table2q, quantile(M2Z(:,1),p1), quantile(M2Z(:,1),p2), quantile(M2Z(:,2),p1), quantile(M2Z(:,2),p2), quantile(M2Z(:,3),p1), quantile(M2Z(:,3),p2),];
Table2q = [Table2q, quantile(M2C(:,1),p1), quantile(M2C(:,1),p2), quantile(M2C(:,2),p1), quantile(M2C(:,2),p2), quantile(M2C(:,3),p1), quantile(M2C(:,3),p2),];
Table2q = [Table2q, quantile(M2W(:,1),p1), quantile(M2W(:,1),p2), quantile(M2W(:,2),p1), quantile(M2W(:,2),p2), quantile(M2W(:,3),p1), quantile(M2W(:,3),p2),];
Table3q = [];
Table3q = [Table3q, quantile(M3Z(:,1),p1), quantile(M3Z(:,1),p2), quantile(M3Z(:,2),p1), quantile(M3Z(:,2),p2), quantile(M3Z(:,3),p1), quantile(M3Z(:,3),p2),];
Table3q = [Table3q, quantile(M3C(:,1),p1), quantile(M3C(:,1),p2), quantile(M3C(:,2),p1), quantile(M3C(:,2),p2), quantile(M3C(:,3),p1), quantile(M3C(:,3),p2),];
Table3q = [Table3q, quantile(M3W(:,1),p1), quantile(M3W(:,1),p2), quantile(M3W(:,2),p1), quantile(M3W(:,2),p2), quantile(M3W(:,3),p1), quantile(M3W(:,3),p2),];

disp('0 vs 1 ME')
set = 1
disp(['Zhang, set =', int2str(set)',': ' ,num2str(signrank(M1Z(:,1),M1Z(:,2)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(signrank(M1C(:,1),M1C(:,2)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(signrank(M1W(:,1),M1W(:,2)))])
set = 2
disp(['Zhang, set =', int2str(set)',': ' ,num2str(signrank(M2Z(:,1),M2Z(:,2)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(signrank(M2C(:,1),M2C(:,2)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(signrank(M2W(:,1),M2W(:,2)))])
set = 3
disp(['Zhang, set =', int2str(set)',': ' ,num2str(signrank(M3Z(:,1),M3Z(:,2)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(signrank(M3C(:,1),M3C(:,2)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(signrank(M3W(:,1),M3W(:,2)))])
%{
disp('0 vs 1 SD')
set = 1
disp(['Zhang, set =', int2str(set)',': ' ,num2str(signrank(S1Z(:,1),S1Z(:,2)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(signrank(S1C(:,1),S1C(:,2)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(signrank(S1W(:,1),S1W(:,2)))])
set = 2
disp(['Zhang, set =', int2str(set)',': ' ,num2str(signrank(S2Z(:,1),S2Z(:,2)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(signrank(S2C(:,1),S2C(:,2)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(signrank(S2W(:,1),S2W(:,2)))])
set = 3
disp(['Zhang, set =', int2str(set)',': ' ,num2str(signrank(S3Z(:,1),S3Z(:,2)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(signrank(S3C(:,1),S3C(:,2)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(signrank(S3W(:,1),S3W(:,2)))])
%}
disp('0 vs 0.85 ME')
set = 1
disp(['Zhang, set =', int2str(set)',': ' ,num2str(signrank(M1Z(:,1),M1Z(:,3)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(signrank(M1C(:,1),M1C(:,3)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(signrank(M1W(:,1),M1W(:,3)))])
set = 2
disp(['Zhang, set =', int2str(set)',': ' ,num2str(signrank(M2Z(:,1),M2Z(:,3)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(signrank(M2C(:,1),M2C(:,3)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(signrank(M2W(:,1),M2W(:,3)))])
set = 3
disp(['Zhang, set =', int2str(set)',': ' ,num2str(signrank(M3Z(:,1),M3Z(:,3)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(signrank(M3C(:,1),M3C(:,3)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(signrank(M3W(:,1),M3W(:,3)))])
%{
disp('0 vs 0.85 SD')
set = 1
disp(['Zhang, set =', int2str(set)',': ' ,num2str(signrank(S1Z(:,1),S1Z(:,3)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(signrank(S1C(:,1),S1C(:,3)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(signrank(S1W(:,1),S1W(:,3)))])
set = 2
disp(['Zhang, set =', int2str(set)',': ' ,num2str(signrank(S2Z(:,1),S2Z(:,3)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(signrank(S2C(:,1),S2C(:,3)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(signrank(S2W(:,1),S2W(:,3)))])
set = 3
disp(['Zhang, set =', int2str(set)',': ' ,num2str(signrank(S3Z(:,1),S3Z(:,3)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(signrank(S3C(:,1),S3C(:,3)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(signrank(S3W(:,1),S3W(:,3)))])
%}



disp('1 vs 0.85')
set = 1
disp(['Zhang, set =', int2str(set)',': ' ,num2str(signrank(M1Z(:,2),M1Z(:,3)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(signrank(M1C(:,2),M1C(:,3)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(signrank(M1W(:,2),M1W(:,3)))])
set = 2
disp(['Zhang, set =', int2str(set)',': ' ,num2str(signrank(M2Z(:,2),M2Z(:,3)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(signrank(M2C(:,2),M2C(:,3)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(signrank(M2W(:,2),M2W(:,3)))])
set = 3
disp(['Zhang, set =', int2str(set)',': ' ,num2str(signrank(M3Z(:,2),M3Z(:,3)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(signrank(M3C(:,2),M3C(:,3)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(signrank(M3W(:,2),M3W(:,3)))])
%{
disp('1 vs 0.85 SD')
set = 1
disp(['Zhang, set =', int2str(set)',': ' ,num2str(signrank(S1Z(:,2),S1Z(:,3)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(signrank(S1C(:,2),S1C(:,3)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(signrank(S1W(:,2),S1W(:,3)))])
set = 2
disp(['Zhang, set =', int2str(set)',': ' ,num2str(signrank(S2Z(:,2),S2Z(:,3)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(signrank(S2C(:,2),S2C(:,3)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(signrank(S2W(:,2),S2W(:,3)))])
set = 3
disp(['Zhang, set =', int2str(set)',': ' ,num2str(signrank(S3Z(:,2),S3Z(:,3)))])
disp(['Carlos, set =' ,int2str(set),': ' ,num2str(signrank(S3C(:,2),S3C(:,3)))])
disp(['Martinez, set =' ,int2str(set),': ',num2str(signrank(S3W(:,2),S3W(:,3)))])

%}