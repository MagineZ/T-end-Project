function  qtvisu2c(dataname)
% Visualize automatically located T-wave ends and annotations by q1c
%
% Author: Qinghua Zhang
% Copyright 2005 INRIA

%Data directory path
QTdir = qtdatapath;

[annodata, Tendindq1c] = DataConvert(dataname, QTdir);

CompareWith1Cardios(annodata, Tendindq1c);
  
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
leftmargin = 1200;  
if dataind1>leftmargin
  dataind1 = dataind1-leftmargin;
else
  leftmargin = dataind1-1;
  dataind1 = 1;
end
dataind2 = find(round(alldata(:,1)*250)==Tendindfin);  
dataind2 = min(dataind2+1200, size(alldata,1));

% Shift Tendindq1c in agreement with annodata
Tendindq1c = Tendindq1c -(Tendinddebut-leftmargin) + 1;

annodata = alldata(dataind1:dataind2,:);    

%================
function  CompareWith1Cardios(annodata, Tendindq1c)
 
fs = 250;
Ts = 1/fs;

nexp = length(Tendindq1c);

Tend0 = twaveend(annodata(:,2), fs);
Tend1 = twaveend(annodata(:,3), fs);

nTend0 = length(Tend0);
nTend1 = length(Tend1);

Tcmp0 = zeros(nexp,1);
Tcmp1 = zeros(nexp,1);

for ks = 1:nexp
  [dum, kc] = min(abs(Tend0-Tendindq1c(ks)));
  Tcmp0(ks) = Tend0(kc);
  [dum, kc] = min(abs(Tend1-Tendindq1c(ks)));
  Tcmp1(ks) = Tend1(kc);  
end

if sum((Tcmp0-Tendindq1c).^2)<sum((Tcmp1-Tendindq1c).^2) 
  bestlead = 0;
else
  bestlead = 1;
end
disp(['Best lead: ' num2str(bestlead)])

clf;
% Plot with drivation 0
subplot(2,1,1)
s0 = annodata(:,2);
  
plotend = max(Tendindq1c) + 1200;
plotend = min(plotend, length(s0));
 
plot((1:plotend)*Ts,s0(1:plotend));
av = axis;
hold on
plot([1;1]*Tend0'*Ts,[ones(1,nTend0)*av(3);ones(1,nTend0)*av(4)],'g');
plot([1;1]*Tendindq1c'*Ts, [ones(1,nexp)*(av(3:4)*[0.5;0.5]); ones(1,nexp)*av(3)],'r')
hold off
  
% Plot with drivation 1
subplot(2,1,2)
s0 = annodata(:,3);
 
plot((1:plotend)*Ts, s0(1:plotend));
av = axis;
hold on
plot([1;1]*Tend1'*Ts,[ones(1,nTend1)*av(3);ones(1,nTend1)*av(4)],'g');
plot([1;1]*Tendindq1c'*Ts, [ones(1,nexp)*(av(3:4)*[0.5;0.5]); ones(1,nexp)*av(3)],'r')
hold off

% Activate context menu for equal xlimits
subzoom init
  
   