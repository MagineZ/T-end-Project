function  qtvisu(dataname)
%QTVISU Visualize automatically located T-wave ends and annotations by q1c & q2c
%qtvisu dataname or qtvisu('dataname') visualize the result for record with
%the file named as dataname.
%
% Author: Qinghua Zhang
% Copyright 2005 INRIA
 
if nargin<1
  error('The data file name should be supplied.')
end

%Data directory path
QTdir = qtdatapath;
if exist([QTdir, dataname, '.Tendq2c'], 'file')
  qtvisu2c(dataname);
elseif exist([QTdir, dataname, '.Tendq1c'], 'file')
  fid = fopen([QTdir, dataname, '.Tendq1c']);
  timecol = textscan(fid,'%s%f%*[^\n]');
  fclose(fid);
  Tendindq1c = timecol{2}; 
  if isempty(Tendindq1c)
    error('No T-wave end annotation for this record.');
  end 
  qtvisu1c(dataname);
else
  errstr = sprintf('The data directory path defined in qtdatapath.m may be incorrect.\nIf it is correct, ');
  error([errstr mfilename ' must be called for record annotated by q1c and/or q2c.'])
end
