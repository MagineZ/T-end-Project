clear; close all ;
cd('/gtmp/peichunsu/Tend data');
table = readtable('MMC_QTauto_QTmanual_allevents.xls');
case_col = {'N733539765','N733549764', 'N733559763','N733569762','N733579761'};
for i = 3%1:5
caseID  = case_col{i};
sss = ismember(table.DeviceSN,caseID);
table1 = table(sss,:);
manual_st = sort(unique(table1.eventDateTime));
seg_col = [];
Tend_col = {};
for j = 1:length(manual_st)
    sss = ismember(table1.eventDateTime,manual_st(j));
    sub_table = table1(sss,:);
    Tends = sub_table.QToffset;
    Tends_manual = sub_table.QTmanoffset;
    Tends(Tends_manual ~= 0) = Tends_manual(Tends_manual ~= 0);
    Tend_col{j} = Tends;
    segID = sub_table.FileNum(1);
    seg_col = [seg_col,segID];
end 
save([caseID,'_ann_1110'],'Tend_col','seg_col');

end