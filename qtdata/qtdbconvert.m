%QTDBCONVERT: convert QT database files to text files
% Read from .dat, .hea, .q1c, .q2c files
% Write to  .txt, .Tendq1c, .Tendq2c, .allq1c, .allq2c
%
% Place this file and the .dat, .hea, .q1c, .q2c files from 
% PhysioNet QT database (available at http://www.physionet.org) 
% in the same directory. Make sure that the commands rdann and rdsamp
% are available (see the WFDB package from PhysioNet). Go to this
% directory in Matlab, and run this file.
% Then the original binary data are converted into text files.
%
% Author: Qinghua Zhang
% Copyright 2005 INRIA


dst=dir('*.q1c');

fichcells = {dst.name};

nfich = length(fichcells);

status = 0;

for kf=1:nfich
  nomk = fichcells{kf}(1:end-4);
  disp(nomk)
  
  [st,rt] = system(['rdann -r ', nomk, ' -a q1c -p '')'' -n 2 >! ', nomk, '.Tendq1c']);
  status = status + st;
  [st,rt] = system(['rdann -r ', nomk, ' -a q1c >! ', nomk, '.allq1c']);
  status = status + st;
  if exist([nomk, '.q2c'], 'file')
    [st,rt] = system(['rdann -r ', nomk, ' -a q2c -p '')'' -n 2 >! ', nomk, '.Tendq2c']);
    status = status + st;
    [st,rt] = system(['rdann -r ', nomk, ' -a q2c >! ', nomk, '.allq2c']);
    status = status + st;
  end
  if status
     error('failed when calling rdann.')
  end
  
  % Determine beginning and end of data conversion
  
  fid = fopen([nomk, '.allq1c']);
  timecol1 = textscan(fid,'%s%*[^\n]');
  timecol1 = timecol1{1};
  fclose(fid);
 
  if exist([nomk, '.allq2c'], 'file')
    fid = fopen([nomk, '.allq2c']);
    timecol2 = textscan(fid,'%s%*[^\n]');
    timecol2 = timecol2{1};
    fclose(fid);
  else
    timecol2 = timecol1;
  end
  
  sec1 = str2num(timecol1{1}(end-5:end));
  minu1 =  str2num(timecol1{1}(1:end-7));
  sec2 = str2num(timecol2{1}(end-5:end));
  minu2 =  str2num(timecol2{1}(1:end-7));
  [dum, ind] = min([minu1*60+sec1 minu2*60+sec2]);
  secsec = [sec1 sec2];
  minuminu = [minu1 minu2];
  sec = secsec(ind);
  minu = minuminu(ind);
  if sec>=6
    sec = sec-6;  
  else
    minu = minu-1;
    sec = sec+54;
  end
  debut = [num2str(minu), ':', num2str(sec)];
  
  sec1 = str2num(timecol1{end}(end-5:end));
  minu1 =  str2num(timecol1{end}(1:end-7));
  sec2 = str2num(timecol2{end}(end-5:end));
  minu2 =  str2num(timecol2{end}(1:end-7));
  [dum, ind] = max([minu1*60+sec1 minu2*60+sec2]);
  secsec = [sec1 sec2];
  minuminu = [minu1 minu2];
  sec = secsec(ind);
  minu = minuminu(ind);
  if sec<54
    sec = sec+6;  
  else
    minu = minu+1;
    sec = sec-54;
  end
  fin = [num2str(minu), ':', num2str(sec)];
  
  % data conversion
  [st,rt] = system(['rdsamp -r ', nomk, ' -p  -f ', debut, ' -t ', fin, ' >! ', nomk, '.txt']);  
  if st
     error('failed when calling rdsamp.')
  end
end

% END of the file
