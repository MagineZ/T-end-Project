%QTDEMOS: demostrations with the PhysioNet QT database
%
% The file qtdatapath.m should be first edited to indicate the directory
% where the data files are located.
%
% Use  qtvisu.m to call T-wave end location algorithm and visualize the
% result.
%
% Author: Qinghua Zhang
% Copyright 2005 INRIA


% Example of ECG record annotated by two cardiologists
qtvisu sel123
% After having zoomed a subplot to examine details, click with the RIGHT
% mouse button on the figure frame to equalize the x-axis of the two
% subplots.


% Example of ECG record annotated by one cardiologist
qtvisu sel47


% Process all the records annotated by two cardiologists by choosing the best lead
% per record.
blprecord2

% Process all the records annotated by two cardiologists by choosing the best lead
% per sample.
blpsample2

% Process all the records annotated by one cardiologist by choosing the best lead
% per record.
blprecord1

% Process all the records annotated by one cardiologist by choosing the best lead
% per sample.
blpsample1


