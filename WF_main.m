
clear; close all ;

%DBfolder = [pwd,'/nifeadb'];
%{
cd('/gtmp/peichunsu/Code Archive3');
addpath ([pwd,'/peak-detector-master'])
addpath([pwd,'/peak-detector-master/sources'])
addpath ([pwd,'/mcode'])
addpath ([pwd,'/decomposition code'])
addpath()
%}
addpath(genpath('/gtmp/peichunsu/Code Archive3'));
addpath(genpath('/gtmp/peichunsu/ecg-kit-0.1.6'));
addpath(genpath('/gtmp/peichunsu/MFEToolbox'))

case_col = {'N733539765','N733549764', 'N733559763','N733569762','N733579761'};
filename_col = {'180510_N733539765','180217_N733549764','180209_N733559763','180319_N733569762', '180220_N733579761'};
cd('/gtmp/peichunsu/Tend data');
table_glucose = readtable('MMC_glucose05022019_a.xls');
table_manual = readtable('MMC_QTauto_QTmanual_allevents.xls');
OptimalShrinkageOpt = 'op'; % Optimal shrinkage norm options: 'fro', 'nuc', 'op'

opt = struct(...
     'JQRS_THRESH', 0.3,... % energy threshold for defining peaks
    'JQRS_REFRAC', 0.25,... % refractory period in seconds
    'JQRS_INTWIN_SZ', 7,...
    'JQRS_WINDOW', 15);


basicTF.win = 500;
basicTF.hop = 20;
basicTF.fs = 100;
basicTF.fr = 0.02;
basicTF.feat = 'SST11';
advTF.num_tap = 1;%num_tap(cc);
advTF.win_type = 'Flattop'; %{'Gaussian','Thomson','multipeak','SWCE'}; %}
%advTF.win_type = 'Hamming'; %{'Gaussian','Thomson','multipeak','SWCE'}; %}
advTF.Smo = 1;
advTF.Rej = 0;
advTF.ths = 1E-6;
advTF.HighFreq = 10/100;
advTF.LowFreq = 0.5/100;
advTF.lpc = 0;
cepR.g = 0.3; % g(cc);%0.06;
cepR.Tc=0;
P.num_s = 1;
P.num_c = 1;
num_tap = 2 ;
lam_curve = 10 ;% lambda for curve extraction (suggest: 10)
lam_beat = 50 ;% lambda for beat tracking (suggest: 20)

qual_all = [];
bsqi_R = [];
bsqi_Tend = [];
q=0;
fs0 = 200;
fs = 1000;
bsqi_Tend_pb = [];
tsqi_all = [];
bsqi_all = [];
for ii = 1:5
    
    
    caseID = case_col{ii};
    fileName = filename_col{ii};
    cd('/gtmp/peichunsu/Tend data');
    load([caseID,'_ann_1110.mat']);
    cd(['/gtmp/peichunsu/Tend data/',caseID]);
    ss = length(dir)-4;
    
    sss = ismember(table_glucose.FileName,fileName );
    sub_table_glucose = table_glucose(sss,:);
    sss = ismember(table_manual.DeviceSN,caseID );
    sub_table_manual = table_manual(sss,:);
    qual = sub_table_glucose.QUAL;
    
    [glu_st,id] = sort(sub_table_glucose.EventTime);
    qual = qual(id);
    qual_all = [qual_all;qual];
   
    num_med_l = round(fs*0.8)+1;
    num_med_s = round(fs*0.4)+1;
    num_med = round(fs*0.1)+1;
    
    [b_hp,a_hp] = butter(4, 0.5/(fs/2), 'high') ;
    [b_lp,a_lp] = butter(4, 30/(fs/2), 'low') ;
  
    wo = 60/(fs/2); bw = wo/35;
    [b_notch1,a_notch1] = iirnotch(wo,bw);
    wo = 50/(fs/2); bw = wo/35;
    [b_notch2,a_notch2] = iirnotch(wo,bw);
    
    final_Rpeaks = {};
    final_Tpeaks = {};
    final_Tend = {};
    
    %}
    for i = 1:length(glu_st)
        
        ann = Tend_col{i}*fs/fs0;
        sss = ismember(sub_table_manual.eventDateTime,glu_st(i));
        cur_table = sub_table_manual(sss,:);
        seg_id = unique(cur_table.FileNum);
        
        ECG = load([pwd,'/',fileName,'_', int2str(seg_id),'.mat']);
        ECG = cell2mat(struct2cell(ECG));
        
        aa = max(round(cur_table.QTonset(1)-60*fs0),1);
        bb = min(round(cur_table.QTonset(end)+60*fs0),length(ECG));
       
    
        sig0 = ECG(aa:bb);
        sig0 = resample(sig0,fs,fs0);
        aa = aa*fs/fs0;
        %[~, Rpeaks0] = twaveend(sig0,fs);
   
        
        x0_real = filtfilt(b_hp,a_hp, sig0) ;
        
        x0_real = filtfilt(b_notch1,a_notch1, x0_real) ;
        x0_real = filtfilt(b_notch2,a_notch2, x0_real) ;
        %x0_real = filtfilt(b_lp,a_lp, x0_real) ;
        x0 = ECG_detrend(x0_real, num_med, num_med, 0,0);
        x0 = x0';x0_real = x0_real';
        %x0_real = ECG_detrend(x0_real, num_med_s, num_med_l,1,0);
        %x0_real = ECG_detrend(x0_real, num_med_l, num_med_l,0,0);
        
        %{
        [~, Rpeaks0] = twaveend(x0',fs);
        
        if median(x0(Rpeaks0)>0)
            [~,Rpeaks] = findrpk_elgendi(x0,fs);
        else
            [~,Rpeaks] = findrpk_elgendi(-x0,fs);
        end
        
        %}
        [~,Rpeaksp] = findrpk_elgendi(x0,fs);
        [~,Rpeaksn] = findrpk_elgendi(-x0,fs);
        if median(x0(Rpeaksp))>median(x0(Rpeaksn))
            Rpeaks = Rpeaksp;
        else
            Rpeaks = Rpeaksn;
        end
        Rpeaks0 = Rpeaks;
        
        %[~, ~,Rpeaks] = extract_Rpeaks_Su(x0,x0_real,fs,basicTF, advTF, cepR, P, lam_curve, lam_beat, 0);
        Rpeaks  = Rpeaks';
        Rpeaks_jqrs =run_qrsdet_by_seg_ali(x0_real,fs,opt);
        uuu=x0_real;
        [Tend1, ~] = twaveend1022(uuu',Rpeaks0,fs);
        
        heasig.freq = fs;
        heasig.nsig = 1;
        heasig.nsamp = length(uuu);
        position0 = wavedet_3D_pcs(uuu',Rpeaks0,fs);
        mT0 = median(position0.Toff(~isnan(position0.Toff))-position0.qrs(~isnan(position0.Toff)));
        position0.Toff(isnan(position0.Toff)) = position0.qrs(isnan(position0.Toff))+floor(mT0);
        Tend2 = position0.Toff;
        sqi_T = [];
        sqi_R = [];
        sqi_wl = 5;
        bsqi_wl=10;
        
        for j = 1:length(Rpeaks)
            Tend_x = Tend1(Tend1>=Rpeaks(j)-sqi_wl*fs);Tend_x = Tend_x(Tend_x<=Rpeaks(j)+sqi_wl*fs);
            Tend_y = Tend2(Tend2>=Rpeaks(j)-sqi_wl*fs);Tend_y = Tend_y(Tend_y<=Rpeaks(j)+sqi_wl*fs);
            R_x = Rpeaks(Rpeaks>=Rpeaks(j)-bsqi_wl*fs);R_x = R_x(R_x<=Rpeaks(j)+bsqi_wl*fs);
            R_y = Rpeaks_jqrs(Rpeaks_jqrs>=Rpeaks(j)-bsqi_wl*fs);R_y = R_y(R_y<=Rpeaks(j)+bsqi_wl*fs);
            sqi_R = [sqi_R, bsqi(R_x,R_y,0.05,fs)];
            sqi_T =  [sqi_T,bsqi(Tend_x, Tend_y, 0.05,fs)];
        end
        
        bsqi_Tend = [bsqi_Tend,median(sqi_T)];
        %A = x0_real;
        %B = x0;
        %beats_todo = Rpeaks(Rpeaks>50*fs);beats_todo = beats_todo(beats_todo<=length(x0)-fs*50);
        beats_todo = Rpeaks;
        num_nonlocal = min(length(beats_todo),20);
        [x0_real,bsqi_Tend_pb0,Rpeaks,sqi_T,sqi_R] = ECG_shrinkage0_qual( x0_real,Rpeaks',beats_todo,sqi_T,sqi_R,q,num_nonlocal, OptimalShrinkageOpt, 1,0);
        nTend = length(ann);
        for ks = 1:nTend
            [~, kc] = min(abs(Rpeaks-ann(ks)));
            sqi0 = bsqi_Tend_pb0(kc);
            bsqi_Tend_pb = [bsqi_Tend_pb,sqi0];
            tsqi_all = [tsqi_all,sqi_T(kc)];
            bsqi_all = [bsqi_all,sqi_R(kc)];
        end
         
        uuu=x0_real;
        [Tend1, ~] = twaveend1022(uuu',Rpeaks0,fs);
        
        heasig.freq = fs;
        heasig.nsig = 1;
        heasig.nsamp = length(uuu);
        position0 = wavedet_3D_pcs(uuu',Rpeaks0,fs);
        mT0 = median(position0.Toff(~isnan(position0.Toff))-position0.qrs(~isnan(position0.Toff)));
        position0.Toff(isnan(position0.Toff)) = position0.qrs(isnan(position0.Toff))+floor(mT0);
        Tend2 = position0.Toff;
        
        
        mT = median(position0.T(~isnan(position0.T))-position0.qrs(~isnan(position0.T)));
        position0.T(isnan(position0.T)) = position0.qrs(isnan(position0.T))+floor(mT);
        Tpeaks = position0.T;
        if isnan(Tpeaks)
            Tend3 = Tpeaks;
        else
            Tend3 = detect_Tend_Carlos(uuu,Tpeaks,fs);
        end
        final_Rpeaks{i} = Rpeaks+aa-1;
        final_Tpeaks{i} = Tpeaks+aa-1;
        final_Tend{i,1} = Tend1+aa-1;
        final_Tend{i,2} = Tend2+aa-1;
        final_Tend{i,3} = Tend3+aa-1;
        
        
    end
    cd('/gtmp/peichunsu/Tend data');
    
    save([caseID,'_locs_0412_shrinkage_sqi_',int2str(fs),'_' num2str(q*100)], 'final_Rpeaks','final_Tpeaks','final_Tend','qual_all','bsqi_Tend','bsqi_Tend_pb','tsqi_all','bsqi_all');
    
    
end
