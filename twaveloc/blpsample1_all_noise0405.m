function   [mmean1,mstd1, mmean2,mstd2] =  blpsample1_all_0312(snr,noise_mth,Tend_mth,sqi_th)
% Best lead per sample, 105 records annotated by one cardiologist
%
% Author: Qinghua Zhang
% Copyright 2005 INRIA

%addpath(genpath('D:\duke box\fECG_project\Code Archive3\Code Archive3'))
%addpath(genpath('D:\duke box\fECG_project\ecg-kit-0.1.6\ecg-kit-0.1.6\common'));
%addpath(genpath('D:\duke box\fECG_project\MFEToolbox'))
%Data directory path
QTdir = qtdatapath;
OptimalShrinkageOpt = 'nuc'; % Optimal shrinkage norm options: 'fro', 'nuc', 'op'
fs = 250;
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

for kd=1:length(dnames)
    dataname = dnames{kd};
    % disp(['Processing record number ', num2str(kd), ', file name: ', dataname])
    
    if (strcmp(dataname, 'sel35') | strcmp(dataname, 'sel37'))
        meanvec1(kd) = 0;
        stdvec1(kd) = 0;
        ncycvec1(kd) = 0;
        continue
    end
    
    [annodata, Tendindq1c] = DataConvert(dataname, QTdir);
    annodata_sh = annodata;
    sig1 = annodata(:,2); sig2 = annodata(:,3);
    x0_real = filtfilt(b_hp,a_hp, sig1') ;
    %x0_real = filtfilt(b_lp,a_lp, x0_real) ;
    x0_real = filtfilt(b_notch1,a_notch1, x0_real) ;
    x0_real = filtfilt(b_notch2,a_notch2, x0_real) ;
    x0 = ECG_detrend(x0_real,num_med,num_med,0,0);
    
    x1_real = filtfilt(b_hp,a_hp, sig2') ;
    %x1_real = filtfilt(b_lp,a_lp, x1_real) ;
    x1_real = filtfilt(b_notch1,a_notch1, x1_real) ;
    x1_real = filtfilt(b_notch2,a_notch2, x1_real) ;
    x1 = ECG_detrend(x1_real,num_med,num_med,0,0);
    %[~,Rpeaks1] = twaveend(sig1,fs);  [~,Rpeaks2] = twaveend(sig2,fs);
    %{
    [~,Rpeaks01] = twaveend(sig1,fs);  [~,Rpeaks02] = twaveend(sig2,fs);
    if median(sig1(Rpeaks01))>0
        [~,Rpeaks1] = findrpk_elgendi(sig1,fs);
    else
        [~,Rpeaks1] = findrpk_elgendi(-sig1,fs);
    end

    
    if median(sig2(Rpeaks02))>0
        [~,Rpeaks2] = findrpk_elgendi(sig2,fs);
    else
        [~,Rpeaks2] = findrpk_elgendi(-sig2,fs);
    end
    %}
    [~,Rpeaks1p] = findrpk_elgendi(x0,fs);[~,Rpeaks1n] = findrpk_elgendi(-x0,fs);
    Rpeaks1p = Rpeaks1p(Rpeaks1p>0);Rpeaks1n = Rpeaks1n(Rpeaks1n>0);
    if median(x0(Rpeaks1p))>median(x0(Rpeaks1n))
        Rpeaks1 = Rpeaks1p;
    else
        Rpeaks1 =Rpeaks1n;
    end
    Rpeaks01 = Rpeaks1;
    [~,Rpeaks2p] = findrpk_elgendi(x1,fs);[~,Rpeaks2n] = findrpk_elgendi(-x1,fs);
    Rpeaks2p = Rpeaks2p(Rpeaks2p>0);Rpeaks2n = Rpeaks2n(Rpeaks2n>0);
    if median(x1(Rpeaks2p))>median(x1(Rpeaks2n))
        Rpeaks2 = Rpeaks2p;
    else
        Rpeaks2 =Rpeaks2n;
    end
    Rpeaks02 = Rpeaks2;
    %Ramp1 = abs(median(sig1(Rpeaks1))); Ramp2 = abs(median(sig2(Rpeaks2)));
    
    if noise_mth == "arma11"
        noise1 = std(sig1)*ARMA11(length(sig1))/(10^(snr/20)); noise2 = std(sig2)*ARMA11(length(sig2))/(10^(snr/20));
        sig1 = sig1+noise1; sig2 = sig2+noise2;
    elseif noise_mth == "gaussian"
        %sig1 = awgn(annodata(:,2),snr,'measured');
        %sig2 = awgn(annodata(:,3),snr,'measured');
        noise1 = std(sig1)*normrnd(0,1,[1,length(sig1)])/(10^(snr/20)); noise2 = std(sig2)*normrnd(0,1,[1,length(sig2)])/(10^(snr/20));
        sig1 = sig1+noise1'; sig2 = sig2+noise2';
    end
    
    x0_real = filtfilt(b_hp,a_hp, sig1') ;
    %x0_real = filtfilt(b_lp,a_lp, x0_real) ;
    x0_real = filtfilt(b_notch1,a_notch1, x0_real) ;
    x0_real = filtfilt(b_notch2,a_notch2, x0_real) ;
    x0 = ECG_detrend(x0_real,num_med,num_med,0,0);
    %x0_real = ECG_detrend(x0_real,num_med_l,num_med_l,0,0);
    
    
    
    x1_real = filtfilt(b_hp,a_hp, sig2') ;
    %x1_real = filtfilt(b_lp,a_lp, x1_real) ;
    x1_real = filtfilt(b_notch1,a_notch1, x1_real) ;
    x1_real = filtfilt(b_notch2,a_notch2, x1_real) ;
    x1 = ECG_detrend(x1_real,num_med,num_med,0,0);
    %x1_real = ECG_detrend(x1_real,num_med_l,num_med_l,0,0);
    
    annodata(:,2) = x0_real;
    annodata(:,3) = x1_real;
    Tend0_Z = twaveend1022(annodata(:,2), Rpeaks01,fs);
    heasig.freq = fs;
    heasig.nsig =1;
    heasig.nsamp = length(annodata(:,2));
    position0 = wavedet_3D_pcs(annodata(:,2), Rpeaks01,fs );
    %position0 = wavedet_3D(annodata(:,2), [],heasig);
    mT0 = median(position0.Toff(~isnan(position0.Toff))-position0.qrs(~isnan(position0.Toff))) ;
    position0.Toff(isnan(position0.Toff)) = position0.qrs(isnan(position0.Toff)) + floor(mT0);
    Tend0_W = position0.Toff;
    sqi_T0 = [];
    sqi_wl = 5;
    for j = 1:length(Rpeaks1)
        Tend_x = Tend0_Z(Tend0_Z>=Rpeaks1(j)-sqi_wl*fs);Tend_x = Tend_x(Tend_x<=Rpeaks1(j)+sqi_wl*fs);
        Tend_y = Tend0_W(Tend0_W>=Rpeaks1(j)-sqi_wl*fs);Tend_y = Tend_y(Tend_y<=Rpeaks1(j)+sqi_wl*fs);    
        sqi_T0 =  [sqi_T0,bsqi(Tend_x, Tend_y, 0.05,fs)];
    end
  
    beats_todo = Rpeaks1;
    num_nonlocal = min(length(Rpeaks1),20);
    x0_real = ECG_shrinkage0_qual1( x0_real,Rpeaks1,beats_todo,sqi_T0,sqi_th,num_nonlocal, OptimalShrinkageOpt, 1,0);
    annodata_sh(:,2) = x0_real;
    
    
    Tend1_Z = twaveend1022(annodata(:,3),Rpeaks02, fs);
    position1 = wavedet_3D_pcs(annodata(:,3), Rpeaks02,fs );
    mT1 = median(position1.Toff(~isnan(position1.Toff))-position1.qrs(~isnan(position1.Toff))) ;
    position1.Toff(isnan(position1.Toff)) = position1.qrs(isnan(position1.Toff)) + floor(mT1);
    Tend1_W = position1.Toff;
    sqi_T1 = [];
    for j = 1:length(Rpeaks2)
        Tend_x = Tend1_Z(Tend1_Z>=Rpeaks2(j)-sqi_wl*fs);Tend_x = Tend_x(Tend_x<=Rpeaks2(j)+sqi_wl*fs);
        Tend_y = Tend1_W(Tend1_W>=Rpeaks2(j)-sqi_wl*fs);Tend_y = Tend_y(Tend_y<=Rpeaks2(j)+sqi_wl*fs);    
        sqi_T1 =  [sqi_T1,bsqi(Tend_x, Tend_y, 0.05,fs)];
    end
   
    beats_todo = Rpeaks2;
    num_nonlocal = min(length(Rpeaks2),20);
    x1_real = ECG_shrinkage0_qual1( x1_real,Rpeaks2,beats_todo,sqi_T1,sqi_th,num_nonlocal, OptimalShrinkageOpt, 1,0);
    annodata_sh(:,3) = x1_real;
    

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
    
end

uerr1 = cell2mat(errcell1);
umean1 = mean(uerr1)*4;
ustd1 = std(uerr1)*4;

uerr2 = cell2mat(errcell1_sh);
umean2 = mean(uerr2)*4;
ustd2 = std(uerr2)*4;

mmean1 = nanmean(meanvec1)*4;
mstd1 = nanmean(stdvec1)*4;
mmean2 = nanmean(meanvec2)*4;
mstd2 = nanmean(stdvec2)*4;

wmean1 = 4*sum(meanvec1.*ncycvec1)/sum(ncycvec1);
wstd1 = 4*sum(stdvec1.*ncycvec1)/sum(ncycvec1);

ncyc1 = sum(ncycvec1);

wmean2 = 4*sum(meanvec2.*ncycvec2)/sum(ncycvec2);
wstd2 = 4*sum(stdvec2.*ncycvec2)/sum(ncycvec2);

ncyc2 = sum(ncycvec2);

%disp(['Auto to Cardio1,', num2str(ncyc1) , ' beats.', ' Result: ', num2str(mmean1), '+/-', num2str(mstd1), ...
%    '  weighted: ', num2str(wmean1), '+/-', num2str(wstd1)])
%disp(['Auto to Cardio1, unique set: ', num2str(umean1), '+/-', num2str(ustd1)])
disp(['Auto to Cardio1(original),', num2str(ncyc1) , ' beats.', ' Result: ', num2str(mmean1), '+/-', num2str(mstd1)]);
disp(['Auto to Cardio1(shrinkage),', num2str(ncyc2) , ' beats.', ' Result: ', num2str(mmean2), '+/-', num2str(mstd2)]);


%Auto to Cardio1,3542 beats. Result: 0.31214+/-17.4331  weighted: 0.37945+/-18.0463
%Auto to Cardio1, unique set: 0.37945+/-38.3476
%Elapsed time is 296.079619 seconds.

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

%x0 = annodata(:,2);
%plot(x0);
%pause;
%}
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

function err1 = CompareWith2Cardios2(annodata, Tendindq1c,Rpeaks0, Rpeaks1)
% Compare automatically located T-wave ends with Tendindq1c, Tendindq2c

fs = 250; % sampling frequency

% compute T-wave ends
%[~,Rpeaks0] = twaveend(annodata(:,2), fs);
%[~,Rpeaks1] = twaveend(annodata(:,3), fs);
%mat2wfdb(annodata(:,2),dataname,fs);
%wrann(dataname,'ann',Rpeaks0,repmat('R',length(Rpeaks0))');
%ecgpuwave(dataname,'atr2',[],[],'ann',0,1);
%Tpeaks0 = rdann(dataname,'atr2',[],[],[],'t');
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

%plot(annodata(:,2))
%pause
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







