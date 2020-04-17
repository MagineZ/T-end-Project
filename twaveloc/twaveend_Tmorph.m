function [Tends, Tpeaks,Rpeaks, asum] = twaveend_Tmorph(s0, fs, swin, mthrld)
%TWAVEEND:  T-wave end location
% Tend = TWAVEEND(S, SF)
% S: ECG signal
% FS: sampling frequency
% Tend: T-wave end indices
% TWAVEEND(S, SF, SWIN, MTHRLD) with optinonal arguments for
% sliding window width and morphology threshold specification.
% SWIN is given as the number of sample points.
% MTHRLD is either a positive number or a character with 'p' for positive
% T-wave or 'n' for negative T-wave.
%
% Author: Qinghua Zhang
% Copyright 2005 INRIA


fratio = fs/250; % Frequence ratio w.r.t. 250hz.

nin = nargin;

if nin<2
  error('Too few input arguments.')
end

if nin<3 | isempty(swin)
  swin = 32*fratio;
end

if nin<4 | isempty(mthrld)
  mthrld = 6;
end

N=length(s0);

[Rpeaks, s] = rpeak(s0, fs);


ptwin = ceil(4*fratio); 

nT = length(Rpeaks);
if Rpeaks(nT)+200 > N
  nT = nT-1;
end

areavalue = zeros(N,1);
Tends = zeros(nT,1);
asum = zeros(nT,1);
Tpeaks = zeros(nT,1);
%xm = 220;
%ymd=70;
%ymu=145;
%%Left part
%ald=0.15;
%bld=ymd-ald*xm;
%alu=0.7;
%blu=ymu-alu*xm;
%%Right part
%ard=0.0;
%brd=ymd-ard*xm;
%aru=0.20;
%bru=ymu-aru*xm;

% T-wave end search interval parameters
ald = 0.15;
bld = 37*fratio;
alu = 0.7;
blu = -9*fratio;
ard = 0.0;
brd = 70*fratio;
aru = 0.20;
bru = 101*fratio;

for knR=1:nT
  kR=Rpeaks(knR);
  
  if (knR<length(Rpeaks))
    RRk = Rpeaks(knR+1)- kR;
  else
    RRk = 200*fratio;
  end
     
  if RRk<220*fratio
    
    minRtoT = floor(ald*RRk+bld); 
    maxRtoT = ceil(alu*RRk+blu);
      
  else
    minRtoT = floor(ard*RRk+brd); 
    maxRtoT = ceil(aru*RRk+bru);
  end
  
  leftbound = kR+minRtoT;
  rightbound = kR+maxRtoT;
    
  rightbound = min(rightbound, N-ptwin);
  leftbound = min(leftbound, rightbound);
  
  if knR<length(Rpeaks) & rightbound>Rpeaks(knR+1)
     rightbound = Rpeaks(knR+1);
  end
  
  % Compute the area indicator
  for kT=leftbound:rightbound
    %cutlevel = mean(s((kT-ptwin):(kT+ptwin)));
    cutlevel = sum(s((kT-ptwin):(kT+ptwin)),1)/(ptwin*2+1);
    corsig = s((kT-swin+1):kT) -  cutlevel;
    areavalue(kT) = sum(corsig,1);
  end
  Tval = areavalue(leftbound:rightbound);
  
  if isnumeric(mthrld) | mthrld(1)=='p'
    [dum, maxind] = max(Tval);
  end
  if isnumeric(mthrld) | mthrld(1)=='n'
    [duminv, maxindinv] = max(-Tval);
  end
  
  if ischar(mthrld)
    if mthrld(1)=='n'
      maxind = maxindinv;
      dum = duminv;
    end
  else  
    if maxind<maxindinv
      leftind = maxind;
      rightind = maxindinv;
      leftdum = dum;
      rightdum = duminv;
      lefts = 1;
      rights = -1;
    else
      leftind = maxindinv;
      rightind = maxind;
      leftdum = duminv;
      rightdum = dum;
      lefts = -1;
      rights = 1;
    end
  
    if leftdum>mthrld*rightdum 
      maxind = leftind;
      dum = leftdum;
      sign = lefts;
      [~,Tind] = max(sign*(s0(leftbound:(leftbound+maxind-1))));
       Tpeaks(knR) = leftbound + Tind-1;
    
    else
      maxind = rightind;
      dum = rightdum;
      sign = rights;
      [~,Tind] = max(sign*(s0((leftbound+rightind-60):(leftbound+rightind-1))));
      Tpeaks(knR) = leftbound+rightind-60 + Tind-1;
  
    end
  end
 
  Tends(knR) = maxind + leftbound - 1;
  asum(knR) = dum;
  %[~,Tind] = max(sign*(s0(leftbound:(leftbound+maxind-1))));
 
end

