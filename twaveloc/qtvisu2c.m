function  qtvisu2c(dataname)
% Visualize automatically located T-wave ends and annotations by q1c & q2c
%
% Author: Qinghua Zhang
% Copyright 2005 INRIA

%Data directory path
QTdir = qtdatapath;
 
[annodata, Tendindq1c, Tendindq2c] = DataConvert(dataname, QTdir);

CompareWith2Cardios(annodata, Tendindq1c, Tendindq2c);
  
%==========================
function [annodata, Tendindq1c, Tendindq2c] = DataConvert(dataname, QTdir)

alldata = load([QTdir, dataname, '.txt']);

fid = fopen([QTdir, dataname, '.Tendq1c']);
timecol = textscan(fid,'%s%f%*[^\n]');
fclose(fid);
Tendindq1c = timecol{2};  

fid = fopen([QTdir, dataname, '.Tendq2c']);
timecol = textscan(fid,'%s%f%*[^\n]');
fclose(fid);
Tendindq2c = timecol{2};

if isempty(Tendindq2c)
  error([mfilename ' must be called for record annotated by both q1c and q2c.'])
end

Tendinddebut = min(Tendindq1c(1),Tendindq2c(1));
Tendindfin = max(Tendindq1c(end),Tendindq2c(end));

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
Tendindq1c = Tendindq1c -(Tendinddebut-leftmargin) + 1;
Tendindq2c = Tendindq2c -(Tendinddebut-leftmargin) + 1;

annodata = alldata(dataind1:dataind2,:);   

%================
function CompareWith2Cardios(annodata, Tendindq1c, Tendindq2c)
 
fs = 250;
Ts = 1/fs;
 
Tend0 = twaveend(annodata(:,2), fs);
Tend1 = twaveend(annodata(:,3), fs);

nTend0 = length(Tend0);
nTend1 = length(Tend1);
 
nq1c = length(Tendindq1c);
nq2c = length(Tendindq2c);

Tcmp0 = zeros(nq1c,1);
Tcmp1 = zeros(nq1c,1);
for ks = 1:nq1c
  [dum, kc] = min(abs(Tend0-Tendindq1c(ks)));
  Tcmp0(ks) = Tend0(kc);
  [dum, kc] = min(abs(Tend1-Tendindq1c(ks)));
  Tcmp1(ks) = Tend1(kc);  
end
if sum((Tcmp0-Tendindq1c).^2)<sum((Tcmp1-Tendindq1c).^2) 
  bestleadq1 = 0;
else
  bestleadq1 = 1;
end
Tcmp0 = zeros(nq2c,1);
Tcmp1 = zeros(nq2c,1);
for ks = 1:nq2c
  [dum, kc] = min(abs(Tend0-Tendindq2c(ks)));
  Tcmp0(ks) = Tend0(kc);
  [dum, kc] = min(abs(Tend1-Tendindq2c(ks)));
  Tcmp1(ks) = Tend1(kc);  
end
if sum((Tcmp0-Tendindq2c).^2)<sum((Tcmp1-Tendindq2c).^2) 
  bestleadq2 = 0;
else
  bestleadq2 = 1;
end

disp(['Best lead for q1c: ' num2str(bestleadq1), ...
   ',  best lead for q2c: ' num2str(bestleadq2)])

clf;
% Plot with drivation 0
subplot(2,1,1)
s0 = annodata(:,2);
plotend = max([Tendindq1c;Tendindq1c]) + 1200;
plotend = min(plotend, length(s0));
 
plot((1:plotend)*Ts,s0(1:plotend));
av = axis;
hold on
plot([1;1]*Tend0'*Ts,[ones(1,nTend0)*av(3);ones(1,nTend0)*av(4)],'g');
plot([1;1]*Tendindq1c'*Ts, [ones(1,nq1c)*(av(3:4)*[1;2])/3; ones(1,nq1c)*av(4)],'r')
plot([1;1]*Tendindq2c'*Ts, [ones(1,nq2c)*av(3); ones(1,nq2c)*(av(3:4)*[2;1])/3],'r')
hold off
  
% Plot with drivation 1
subplot(2,1,2)
s0 = annodata(:,3);
 
plot((1:plotend)*Ts,s0(1:plotend));
av = axis;
hold on
plot([1;1]*Tend1'*Ts,[ones(1,nTend1)*av(3);ones(1,nTend1)*av(4)],'g');
plot([1;1]*Tendindq1c'*Ts, [ones(1,nq1c)*(av(3:4)*[1;2])/3; ones(1,nq1c)*av(4)],'r')
plot([1;1]*Tendindq2c'*Ts, [ones(1,nq2c)*av(3); ones(1,nq2c)*(av(3:4)*[2;1])/3],'r')
hold off
  
% Activate context menu for equal xlimits
subzoom init 
 
