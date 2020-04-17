%OTHERDEMOS demonstration with arbitrary ECG signal
%
% Author: Qinghua Zhang
% Copyright 2005 INRIA

% Load any ECG signal. In this example a signal of the QT data base is used.
load([qtdatapath 'sel103.txt'])
s = sel103(:,3);

% Compute T-wave end locations 
tends = twaveend(s, 250);
% Visualize the result
locplot(s, 250, tends)

% Compute R-peak locations 
rps = rpeak(s, 250);
% Visualize the result
locplot(s, 250, rps, tends)