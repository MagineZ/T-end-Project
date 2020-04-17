clear; close all;
%DBfolder = [pwd,'/nifeadb'];
cd('/gtmp/peichunsu/Code Archive3');
addpath ([pwd,'/peak-detector-master'])
addpath([pwd,'/peak-detector-master/sources'])
addpath ([pwd,'/mcode'])
addpath ([pwd,'/decomposition code'])
cd('/gtmp/peichunsu/Tend data');

table = readtable('MMC_QTauto_QTmanual_allevents.xls');
case_col = {'N733539765','N733549764', 'N733559763','N733569762','N733579761'};
filename_col = {'180510_N733539765','180217_N733549764','180209_N733559763','180319_N733569762', '180220_N733579761'};
fs = 250;
fs0=200;
total_me = {[],[],[]};
total_qual = [];
    count = 0;
    for i = 1:5
        
        caseID  = case_col{i};
        load([caseID,'_ann_1110.mat']);
        load([caseID,'_locs_0412_shrinkage_sqi_' int2str(fs) '_0.mat']);
        %load([caseID,'_locs_0329.mat']);
        %load([caseID,'_locs_0219.mat']);
        %load([caseID,'_locs.mat']);x0,x0_real] = ECG_shrinkage0(x0,x0_real,Rpeaks',OptimalShrinkageOpt,1,0);
             

        Ttest_all = {};
        ann_all={};
        total_ann_number = 0;
         
        for j = 1:length(Tend_col)
           count = count+1;
            ann = Tend_col{j}*fs/fs0;
            total_ann_number = total_ann_number + length(ann);
          
                for mth = 1:3
                    
                    %Tends = final_Tend{seg_col(j),mth};
                    Tends = final_Tend{j,mth};
                    %Tends = Tends(Tends>(ann(1)-round(fs/20)));
                    %Tends = Tends(Tends<(ann(end)+round(fs/20)));
                    err = CompareWith2Cardios(Tends,ann);
                    %total_me(count,mth) = nanmean(abs(err))*5;
                    %total_sd(count,mth) = nanstd(abs(err))*5;
                    total_me{mth} = [total_me{mth};abs(err)*1000/fs];
                    
                    
                end
               
         total_qual = [total_qual, qual_all(count)*ones(1,length(err))];
              
        end
        
       
        
    end
    
 Z = total_me{1};     C = total_me{2};    W = total_me{3};
 bsqi_Tend_pb1 = bsqi_Tend_pb(total_qual == 1);
 bsqi_Tend_pb3 = bsqi_Tend_pb(total_qual == 3);
 bsqi_Tend_pb5 = bsqi_Tend_pb(total_qual == 5);
 
 
total_me = {[],[],[]};
total_qual = [];
    count = 0;
    for i = 1:5
        
        caseID  = case_col{i};
        load([caseID,'_ann_1110.mat']);
        load([caseID,'_locs_0412_shrinkage_sqi_' int2str(fs) '_90.mat']);
      
        %load([caseID,'_locs_0329.mat']);
        %load([caseID,'_locs_0219.mat']);
        %load([caseID,'_locs.mat']);x0,x0_real] = ECG_shrinkage0(x0,x0_real,Rpeaks',OptimalShrinkageOpt,1,0);
             

        Ttest_all = {};
        ann_all={};
        total_ann_number = 0;
         
        for j = 1:length(Tend_col)
           count = count+1;
            ann = Tend_col{j}*fs/fs0;
            total_ann_number = total_ann_number + length(ann);
          
                for mth = 1:3
                    
                    %Tends = final_Tend{seg_col(j),mth};
                    Tends = final_Tend{j,mth};
                    %Tends = Tends(Tends>(ann(1)-round(fs/20)));
                    %Tends = Tends(Tends<(ann(end)+round(fs/20)));
                    err = CompareWith2Cardios(Tends,ann);
                    %total_me(count,mth) = nanmean(abs(err))*5;
                    %total_sd(count,mth) = nanstd(abs(err))*5;
                    total_me{mth} = [total_me{mth};abs(err)*1000/fs];
                    
                    
                end
               
         total_qual = [total_qual, qual_all(count)*ones(1,length(err))];
              
        end
        
       
        
    end
    
    
 Z = [Z,total_me{1}];     C = [C,total_me{2}];    W = [W,total_me{3}];
 bsqi_Tend_pb1 = bsqi_Tend_pb(total_qual == 1);
 bsqi_Tend_pb3 = bsqi_Tend_pb(total_qual == 3);
 bsqi_Tend_pb5 = bsqi_Tend_pb(total_qual == 5);
 LL =W(total_qual == 1,: );
 LL = LL(bsqi_Tend_pb1 == 1,:);
 median(LL)
 mad(LL)
 quantile(LL,0.025)
  quantile(LL,0.975)
  save(['WF_stat_' int2str(fs)],'Z','C','W','bsqi_Tend_pb1','bsqi_Tend_pb3','bsqi_Tend_pb5','tsqi_all','bsqi_all')