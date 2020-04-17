function   [uerr1,uerr2,bsqi_Tend_pb, tsqi_set1, bsqi_set1] =  blpsample1_all_sqi_0406(Tend_mth,sqi_th,fs)
% Best lead per sample, 105 records annotated by one cardiologist
%Data directory path
opt = struct(...
    ... % jqrs parameters - the custom peak detector implemented herein
    'JQRS_THRESH', 0.3,... % energy threshold for defining peaks
    'JQRS_REFRAC', 0.25,... % refractory period in seconds
    'JQRS_INTWIN_SZ', 7,...
    'JQRS_WINDOW', 15);

QTdir = qtdatapath;
OptimalShrinkageOpt = 'op'; % Optimal shrinkage norm options: 'fro', 'nuc', 'op'
fs0 = 250;
num_med = round(fs*0.1)+1;
num_med_s =  round(fs*0.4)+1;
num_med_l =  round(fs*0.8)+1;
[b_hp,a_hp] = butter(4, 0.5/(fs/2), 'high') ;
[b_lp,a_lp] = butter(4, 30/(fs/2), 'low') ;
wo = 60/(fs/2);  bw = wo/35;
[b_notch1,a_notch1] = iirnotch(wo,bw);
wo = 50/(fs/2);  bw = wo/35;
[b_notch2,a_notch2] = iirnotch(wo,bw);

% Get the list of data file names
dst = dir([QTdir, '*.Tendq1c']);
fichcells = {dst.name};
nfich = length(fichcells);
dnames = cell(nfich,1);
for kd=1:nfich
    dnames{kd} = fichcells{kd}(1:end-8);
end

errcell1 = cell(nfich,1);
errcell1_sh = cell(nfich,1);
% Memory initialization
meanvec1 = zeros(nfich,1);
stdvec1 = zeros(nfich,1);

ncycvec1 = zeros(nfich,1);
siglength = zeros(nfich,1);
bsqi_Tend_pb = [];
tsqi_ch1 = [];
tsqi_ch2 = [];
bsqi_ch1 = [];
bsqi_ch2 = [];
    
for kd= 1:length(dnames) %31
    dataname = dnames{kd};
    %disp(['Processing record number ', num2str(kd), ', file name: ', dataname])
    
    if (strcmp(dataname, 'sel35') | strcmp(dataname, 'sel37'))   
        meanvec1(kd) = 0;
        stdvec1(kd) = 0;
        ncycvec1(kd) = 0;
        continue
    end
    
    [annodata, Tendindq1c] = DataConvert(dataname, QTdir);
    annodata = [resample(annodata(:,2),fs,fs0),resample(annodata(:,3),fs,fs0)];
    annodata = [ones(length(annodata(:,2)),1),annodata];
    Tendindq1c =  Tendindq1c*fs/fs0;
    annodata_sh = annodata;

    sig1 = annodata(:,2); sig2 = annodata(:,3);
    x1_real = sig1-mean(sig1); x2_real = sig2-mean(sig2);
    x1_real = filtfilt(b_hp,a_hp, x1_real') ;
    x1_real = filtfilt(b_notch1,a_notch1, x1_real) ;
    x1_real = filtfilt(b_notch2,a_notch2, x1_real) ;
    x1 = ECG_detrend(x1_real,num_med,num_med,0,0);
    
   
    x2_real = filtfilt(b_hp,a_hp, x2_real') ;
    x2_real = filtfilt(b_notch1,a_notch1, x2_real) ;
    x2_real = filtfilt(b_notch2,a_notch2, x2_real) ;
    x2 = ECG_detrend(x2_real,num_med,num_med,0,0);
    
    
    [~,Rpeaks1p] = findrpk_elgendi(x1,fs);[~,Rpeaks1n] = findrpk_elgendi(-x1,fs);
    Rpeaks1p = Rpeaks1p(Rpeaks1p>0);Rpeaks1n = Rpeaks1n(Rpeaks1n>0);
    if median(x1(Rpeaks1p))>median(x1(Rpeaks1n))
        Rpeaks1 = Rpeaks1p;
    else
        Rpeaks1 =Rpeaks1n;
    end
    Rpeaks1_jqrs =run_qrsdet_by_seg_ali(x1*10^4,fs,opt);
    

    [~,Rpeaks2p] = findrpk_elgendi(x2,fs);[~,Rpeaks2n] = findrpk_elgendi(-x2,fs);
    Rpeaks2p = Rpeaks2p(Rpeaks2p>0);Rpeaks2n = Rpeaks2n(Rpeaks2n>0);
    if median(x2(Rpeaks2p))>median(x2(Rpeaks2n))
        Rpeaks2 = Rpeaks2p;
    else
        Rpeaks2 =Rpeaks2n;
    end
    Rpeaks2_jqrs =run_qrsdet_by_seg_ali(x2*10^4,fs,opt);
   
    annodata(:,2) = x1_real;
    annodata(:,3) = x2_real;
    Tend0_Z = twaveend1022(annodata(:,2), Rpeaks1,fs);
    heasig.freq = fs;
    heasig.nsig =1;
    heasig.nsamp = length(annodata(:,2));
    position0 = wavedet_3D_pcs(annodata(:,2), Rpeaks1,fs );
    mT0 = median(position0.Toff(~isnan(position0.Toff))-position0.qrs(~isnan(position0.Toff))) ;
    position0.Toff(isnan(position0.Toff)) = position0.qrs(isnan(position0.Toff)) + floor(mT0);
    Tend0_W = position0.Toff;
    sqi_T0 = [];
    sqi_R0 = [];
    sqi_wl = 5;
    bsqi_wl = 10;
    for j = 1:length(Rpeaks1)
        Tend_x = Tend0_Z(Tend0_Z>=Rpeaks1(j)-sqi_wl*fs);Tend_x = Tend_x(Tend_x<=Rpeaks1(j)+sqi_wl*fs);
        Tend_y = Tend0_W(Tend0_W>=Rpeaks1(j)-sqi_wl*fs);Tend_y = Tend_y(Tend_y<=Rpeaks1(j)+sqi_wl*fs);   
        R_x = Rpeaks1(Rpeaks1>=Rpeaks1(j)-bsqi_wl*fs);R_x = R_x(R_x<=Rpeaks1(j)+bsqi_wl*fs);
        R_y = Rpeaks1_jqrs(Rpeaks1_jqrs>=Rpeaks1(j)-bsqi_wl*fs);R_y = R_y(R_y<=Rpeaks1(j)+bsqi_wl*fs);
        sqi_T0 =  [sqi_T0,bsqi(Tend_x, Tend_y, 0.05,fs)];
        sqi_R0 = [sqi_R0, bsqi(R_x,R_y,0.05,fs)];
    end
    
    beats_todo = Rpeaks1;
    num_nonlocal = min(length(Rpeaks1),20);
    [x1_real,bsqi_Tend_pb0,Rpeaks1,sqi_T0,sqi_R0] = ECG_shrinkage0_qual( x1_real,Rpeaks1,beats_todo,sqi_T0,sqi_R0,sqi_th,num_nonlocal, OptimalShrinkageOpt, 1,0);
    annodata_sh(:,2) = x1_real;
    
   
    Tend1_Z = twaveend1022(annodata(:,3),Rpeaks2, fs);
    position1 = wavedet_3D_pcs(annodata(:,3), Rpeaks2,fs );
    mT1 = median(position1.Toff(~isnan(position1.Toff))-position1.qrs(~isnan(position1.Toff))) ;
    position1.Toff(isnan(position1.Toff)) = position1.qrs(isnan(position1.Toff)) + floor(mT1);
    Tend1_W = position1.Toff;
    sqi_T1 = [];
    sqi_R1 = [];
    for j = 1:length(Rpeaks2)
        Tend_x = Tend1_Z(Tend1_Z>=Rpeaks2(j)-sqi_wl*fs);Tend_x = Tend_x(Tend_x<=Rpeaks2(j)+sqi_wl*fs);
        Tend_y = Tend1_W(Tend1_W>=Rpeaks2(j)-sqi_wl*fs);Tend_y = Tend_y(Tend_y<=Rpeaks2(j)+sqi_wl*fs);    
        R_x = Rpeaks2(Rpeaks2>=Rpeaks2(j)-bsqi_wl*fs);R_x = R_x(R_x<=Rpeaks2(j)+bsqi_wl*fs);
        R_y = Rpeaks2_jqrs(Rpeaks2_jqrs>=Rpeaks2(j)-bsqi_wl*fs);R_y = R_y(R_y<=Rpeaks2(j)+bsqi_wl*fs);
        sqi_T1 =  [sqi_T1,bsqi(Tend_x, Tend_y, 0.05,fs)];
        sqi_R1 = [sqi_R1, bsqi(R_x,R_y,0.05,fs)];
    end
    
    beats_todo = Rpeaks2;
    num_nonlocal = min(length(Rpeaks2),20);
    [x2_real,bsqi_Tend_pb1,Rpeaks2,sqi_T1,sqi_R1] = ECG_shrinkage0_qual( x2_real,Rpeaks2,beats_todo,sqi_T1,sqi_R1,sqi_th,num_nonlocal, OptimalShrinkageOpt, 1,0);
    annodata_sh(:,3) = x2_real;
    
   
    nTend = length(Tendindq1c);
    for ks = 1:nTend
        [~, kc] = min(abs(Rpeaks1-Tendindq1c(ks)));
        s0 = bsqi_Tend_pb0(kc);
        tsqi_ch1 = [tsqi_ch1, sqi_T0(kc)];
        bsqi_ch1 = [bsqi_ch1,sqi_R0(kc)];
        [~, kc] = min(abs(Rpeaks2-Tendindq1c(ks)));
        s1 = bsqi_Tend_pb1(kc);
        tsqi_ch2 = [tsqi_ch2, sqi_T1(kc)];
        bsqi_ch2 = [bsqi_ch2,sqi_R1(kc)];
        bsqi_Tend_pb = [bsqi_Tend_pb,(s0) ||(s1)];
    end
    
   
    
    if Tend_mth == "Zhang"
        errcell1{kd} = CompareWith2Cardios1(annodata, Tendindq1c,Rpeaks1,Rpeaks2);
        errcell1_sh{kd} = CompareWith2Cardios1(annodata_sh, Tendindq1c,Rpeaks1,Rpeaks2);
    elseif Tend_mth == "Carlos"
        errcell1{kd} = CompareWith2Cardios2(annodata, Tendindq1c,Rpeaks1,Rpeaks2);
        errcell1_sh{kd} = CompareWith2Cardios2(annodata_sh, Tendindq1c,Rpeaks1,Rpeaks2);
    elseif Tend_mth == "Wavedet"
        errcell1{kd} = CompareWith2Cardios3(annodata, Tendindq1c,Rpeaks1,Rpeaks2);
        errcell1_sh{kd} = CompareWith2Cardios3(annodata_sh, Tendindq1c,Rpeaks1,Rpeaks2);
    end
   
    errcell1{kd} = abs(errcell1{kd});
    errcell1_sh{kd} = abs(errcell1_sh{kd});
    
    meanvec1(kd) = nanmean(errcell1{kd});
    stdvec1(kd) = nanstd(errcell1{kd});
    ncycvec1(kd) = length(errcell1{kd});
    
    meanvec2(kd) = nanmean(errcell1_sh{kd});
    stdvec2(kd) = nanstd(errcell1_sh{kd});
    ncycvec2(kd) = length(errcell1_sh{kd});
    siglength(kd) = length(annodata(:,2));
    
end
tsqi_set1 = [tsqi_ch1;tsqi_ch2];
bsqi_set1 = [bsqi_ch1;bsqi_ch2];
uerr1 = cell2mat(errcell1)*1000/fs;
uerr2 = cell2mat(errcell1_sh)*1000/fs;


%========================================
function [annodata, Tendindq1c] = DataConvert(dataname, QTdir)
% data convertion
% annodata: annotated portion of data
% Tendindq1c: T-wave end indices annotated by q1c

alldata = load([QTdir, dataname, '.txt']);
fid = fopen([QTdir, dataname, '.Tendq1c']);
timecol = textscan(fid,'%s%f%*[^\n]');
fclose(fid);
Tendindq1c = timecol{2};

% Determine beginning and end of annotated data
Tendinddebut =  Tendindq1c(1);
Tendindfin = Tendindq1c(end);
dataind1 = find(round(alldata(:,1)*250)==Tendinddebut);
leftmargin = 60*250;%1200;
if dataind1>leftmargin
    dataind1 = dataind1-leftmargin;
else
    leftmargin = dataind1-1;
    dataind1 = 1;
end
dataind2 = find(round(alldata(:,1)*250)==Tendindfin);
dataind2 = min(dataind2+leftmargin, size(alldata,1));

% Shift Tendindq1c in agreement with annodata
Tendindq1c = Tendindq1c -(Tendinddebut-leftmargin) + 1;

annodata = alldata(dataind1:dataind2,:);

%================
function err1 = CompareWith2Cardios1(annodata, Tendindq1c,Rpeaks1,Rpeaks2)
% Compare automatically located T-wave ends with Tendindq1c, Tendindq2c

fs = 250; % sampling frequency

% compute T-wave ends

%Tend0 = twaveend(annodata(:,2), fs);
%Tend1 = twaveend(annodata(:,3), fs);
Tend0 = twaveend1022(annodata(:,2), Rpeaks1,fs);
Tend1 = twaveend1022(annodata(:,3),Rpeaks2, fs);
%[Tend0 ,T0]= twaveend_Tmorph(annodata(:,2), fs);
%[Tend1 ,T1]= twaveend_Tmorph(annodata(:,3), fs);


%compare with Tendindq1c
nTend = length(Tendindq1c);
Tcmp1 = zeros(nTend,1);
for ks = 1:nTend
    [dum, kc] = min(abs(Tend0-Tendindq1c(ks)));
    Tq1d0 = Tend0(kc);
    [dum, kc] = min(abs(Tend1-Tendindq1c(ks)));
    Tq1d1 = Tend1(kc);
    if abs(Tq1d0-Tendindq1c(ks))<abs(Tq1d1-Tendindq1c(ks))
        Tcmp1(ks) = Tq1d0;
    else
        Tcmp1(ks) = Tq1d1;
    end
    %{
    er = abs(Tcmp1(ks) - Tendindq1c(ks));
    if er >40 && ifplot
        t = length(annodata(:,2))
        plot(annodata(:,2)); hold on;
        plot(annodata(:,3)-.5);
        xline(Tendindq1c(ks),'k');
        xline(Tcmp1(ks),'r');
        hold off;
        xlim([Tcmp1(ks)-250,Tcmp1(ks)+80])
        pause;
        ss = randi(20);
        saveas(gcf,['QT_test_figure_', int2str(ss)],'jpg');
    end
    %}
end

err1 = Tcmp1-Tendindq1c;

function err1 = CompareWith2Cardios2(annodata, Tendindq1c,Rpeaks0, Rpeaks1)
% Compare automatically located T-wave ends with Tendindq1c, Tendindq2c

fs = 250; % sampling frequency

% compute T-wave ends
heasig.freq = fs;
heasig.nsig =1;
heasig.nsamp = length(annodata(:,2));
position0 = wavedet_3D_pcs(annodata(:,2), Rpeaks0,fs );
%position0 = wavedet_3D(annodata(:,2), [],heasig);
mT0 = median(position0.T(~isnan(position0.T))-position0.qrs(~isnan(position0.T))) ;
position0.T(isnan(position0.T)) = position0.qrs(isnan(position0.T)) + floor(mT0);
Tpeaks0 = position0.T;
if isnan(Tpeaks0)
    Tend0 = Tpeaks0;
else
    Tend0 = detect_Tend_Carlos(annodata(:,2)',Tpeaks0,fs);
end
%[~,~,Tpeaks0] = elgendi_pt(annodata(:,2),Rpeaks0,round(55/4),round(110/4));
position1 = wavedet_3D_pcs(annodata(:,3), Rpeaks1,fs );
%position1 = wavedet_3D(annodata(:,3), [],heasig);
mT1 = median(position1.T(~isnan(position1.T))-position1.qrs(~isnan(position1.T))) ;
position1.T(isnan(position1.T)) = position1.qrs(isnan(position1.T)) + floor(mT1);
Tpeaks1 = position1.T;
%[~,Rpeaks1] = twaveend(annodata(:,3), fs);
%[~,~,Tpeaks1] = elgendi_pt(annodata(:,3),Rpeaks1,round(55/4),round(110/4));
if isnan(Tpeaks1)
    Tend1 = Tpeaks1;
else
    Tend1 = detect_Tend_Carlos(annodata(:,3)',Tpeaks1,fs);
end

%compare with Tendindq1c
nTend = length(Tendindq1c);
Tcmp1 = zeros(nTend,1);
for ks = 1:nTend
    [dum, kc] = min(abs(Tend0-Tendindq1c(ks)));
    Tq1d0 = Tend0(kc);
    [dum, kc] = min(abs(Tend1-Tendindq1c(ks)));
    Tq1d1 = Tend1(kc);
    if abs(Tq1d0-Tendindq1c(ks))<abs(Tq1d1-Tendindq1c(ks))
        Tcmp1(ks) = Tq1d0;
    else
        Tcmp1(ks) = Tq1d1;
    end
    
end

err1 = Tcmp1-Tendindq1c;

function err1 = CompareWith2Cardios3(annodata, Tendindq1c,Rpeaks0,Rpeaks1)
% Compare automatically located T-wave ends with Tendindq1c, Tendindq2c

fs = 250; % sampling frequency

% compute T-wave ends

%[~,Rpeaks0] = twaveend(annodata(:,2), fs);
%[~,Rpeaks1] = twaveend(annodata(:,3), fs);
%position0 = wavedet_3D_pcs(annodata(:,2), Rpeaks0,fs );
heasig.freq = fs;
heasig.nsig =1;
heasig.nsamp = length(annodata(:,2));

position0 = wavedet_3D_pcs(annodata(:,2), Rpeaks0,fs );
%position0 = wavedet_3D(annodata(:,2), [],heasig);
mT0 = median(position0.Toff(~isnan(position0.Toff))-position0.qrs(~isnan(position0.Toff))) ;
position0.Toff(isnan(position0.Toff)) = position0.qrs(isnan(position0.Toff)) + floor(mT0);
Tend0 = position0.Toff;
position1 = wavedet_3D_pcs(annodata(:,3), Rpeaks1,fs );
%position1 = wavedet_3D(annodata(:,3), [],heasig);
mT1 = median(position1.Toff(~isnan(position1.Toff))-position1.qrs(~isnan(position1.Toff))) ;
position1.Toff(isnan(position1.Toff)) = position1.qrs(isnan(position1.Toff)) + floor(mT1);
Tend1 = position1.Toff;
%compare with Tendindq1c


nTend = length(Tendindq1c);
Tcmp1 = zeros(nTend,1);
for ks = 1:nTend
    [dum, kc] = min(abs(Tend0-Tendindq1c(ks)));
    Tq1d0 = Tend0(kc);
    [dum, kc] = min(abs(Tend1-Tendindq1c(ks)));
    Tq1d1 = Tend1(kc);
    if abs(Tq1d0-Tendindq1c(ks))<abs(Tq1d1-Tendindq1c(ks))
        Tcmp1(ks) = Tq1d0;
    else
        Tcmp1(ks) = Tq1d1;
    end
    
end

err1 = Tcmp1-Tendindq1c;







