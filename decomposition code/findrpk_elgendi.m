function [rri,locsf,tlocsf,ampsf] = findrpk_elgendi(in,fs,varargin)

% R peak detection based on Elgendi algo
% Name:
%   findrpk_elgendi.m
%
% Description:
%   Extraction of the R-R interval time series based on Elgendi's algo
%
% Input:
%   in:             Input signal, 2-column format
%                    Column 1:   Time stamps
%                    Column 2:   ECG value
%   varargin:       optional display argument
%
% Output:
%   locsf:          indices where the R-peaks are located
%   tlocsf:         corresponding timestamps
%   ampsf:          corresponding amplitudes
%   rri:            R-R interval time series
%
% Toolbox and External Dependencies:
%
% References:
%  M. Elgendi, Fast QRS Detection with an Optimized Knowledge-Based 
%   Method: Evaluation on 11 Standard ECG Databases, PLOS one 8(9):e73557 2013
%
% Author: Dynamical Analysis Laboratory, Ottawa Hospital Research Institute
% Date: June 2015
%


if ~isempty(varargin)
    if size(varargin,2)==1
        if ~isempty(varargin{1})
            opt=varargin{1};
        else
            opt.Elgendi_beat_window_length=ceil(220*fs/360);
            opt.Elgendi_qrs_window_length=ceil(35*fs/360);
            opt.Elgendi_beta=1/8;
            
        end
        displayflag=0;
        detection=[];
    elseif size(varargin,2)==2
        if ~isempty(varargin{1})
            opt=varargin{1};
        else
            opt.Elgendi_beat_window_length=ceil(220*fs/360);
            opt.Elgendi_qrs_window_length=ceil(35*fs/360);
            opt.Elgendi_beta=1/8;
            
        end
        displayflag=varargin{2};
        detection=[];
    else
        if ~isempty(varargin{1})
            opt=varargin{1};
        else
            opt.Elgendi_beat_window_length=ceil(220*fs/360);
            opt.Elgendi_qrs_window_length=ceil(35*fs/360);
            opt.Elgendi_beta=1/8;
            
        end
        displayflag=varargin{2};
        detection=varargin{3};
    end
else
    opt.Elgendi_beat_window_length=ceil(220*fs/360);
    opt.Elgendi_qrs_window_length=ceil(35*fs/360);
    opt.Elgendi_beta=1/8;
    
    displayflag=0;
    detection=[];
end

t1= [1:length(in)] / fs ;
x1 = in ;
wr=[2*8/fs 2*20/fs];
[br,ar]=butter(3,wr);
y=filtfilt(br,ar,x1);
W2=opt.Elgendi_beat_window_length;
W1=opt.Elgendi_qrs_window_length;
beta=opt.Elgendi_beta;
y=y.^2;
MAqrs=movmean(y,W1);
MAbeat=movmean(y,W2);

z=mean(y);
alpha=beta*z;
THR1=MAbeat+alpha;
BOF(length(MAqrs),1)=0;

BOF(MAqrs>THR1)=0.1;

[M,V] = regexp(sprintf('%i',[0 BOF'==0.1]),'1+','match');
Nn=cellfun('length',M);
st=V-1;
et=min([st+Nn; length(x1)*ones(1,length(st))]);
nblocks=numel(st);
locsf=[];
THR2=W1;
for j=1:nblocks
    if Nn(j)>=THR2
        [~,R_peak]=max(x1(st(j):et(j)));
        locsf=[locsf R_peak+st(j)-1];
    end
end

%%
%Additional step NOT in original paper
%Filter peak locations if too close based on refractory period of 250ms
refrac_time=floor(0.25*fs); %
rd=diff(locsf);

indmoins=find(rd<refrac_time);
dis=[];
slopdex=10;
for kkk=1:length(indmoins)
    
    preslope1=in(locsf(indmoins(kkk)))-in(max(1,locsf(indmoins(kkk))-slopdex));
    postslope1=-(in(locsf(indmoins(kkk))+slopdex)-in(locsf(indmoins(kkk))));
    maxslope1=min([preslope1 postslope1]);
    preslopep1=in(locsf(indmoins(kkk)+1))-in(locsf(indmoins(kkk)+1)-slopdex);
    postslopep1=-(in(locsf(indmoins(kkk)+1)+slopdex)-in(locsf(indmoins(kkk)+1)));
    maxslopep1=min([preslopep1 postslopep1]);
    if maxslopep1>maxslope1
        dis=[dis; locsf(indmoins(kkk))];
    else
        dis=[dis; locsf(indmoins(kkk)+1)];
    end
end
dis=unique(dis);
locsf=setdiff(locsf,dis);

%Collect results and outputs
if numel(locsf)>1
    tlocsf=t1(locsf);
    ampsf=x1(locsf);
    tr = t1(locsf);
    rri(:,1) = tr(2:end);
    C = tr(2:end) - tr(1:end-1);
    D = datevec (C); % Convert HR intervals to Matlab date vector object.
    % Convert datevector object to seconds
    yyyy_ss = D(:,1)*365.25*24*60*60*1;
    mm_ss = D(:,2)*30*24*60*60*1;
    dd_ss = D(:,3)*24*60*60*1;
    HH_ss = D(:,4)*60*60*1;
    MM_ss = D(:,5)*60*1;
    SS_ss = D(:,6)*1;
    E =  yyyy_ss + mm_ss + dd_ss + HH_ss + MM_ss + SS_ss;
    rri(:,2) = E;
    
else
    rri=[];
    locsf=[];
    tlocsf=[];
    ampsf=[];
end

SearchLen = round(fs/20);
x_up = in;
for ii = 1:length(locsf)
    [~, idx] = max(x_up(max([locsf(ii)-SearchLen 1]):min([locsf(ii)+SearchLen length(x_up)]))) ;
    locsf(ii) = locsf(ii)-SearchLen+idx-1 ;
   
end

