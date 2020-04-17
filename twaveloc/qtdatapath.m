function QTdir = qtdatapath
%qtdatapath indicate the QT data directory 
%
% The line defining QTdir should be edited by the user to indicate 
% the absolute path of QT data directory.
%
% Author: Qinghua Zhang
% Copyright 2005 INRIA

% Edit the following line to indicate the absolute path of data directory
QTdir ='D:\duke box\fECG_project\ecgtwave\ecgtwave\qtdata\';
% ATTENTION: do not omit the the '/' at the end of the path.

if ~exist('QTdir', 'var')
  error('QTdir is not defined.')
end
if isempty(QTdir) | (QTdir(end)~='/' & QTdir(end)~='\')
  error('QTdir must be a string ended with ''/'' or ''\''.')
end
