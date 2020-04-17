function [PR_DR,QT_DR,ST_DR, FTQRS_DR, FPQRS_DR,FTQRS0,FTQRS,FPQRS0,FPQRS,T0_height,T_height,P0_height,P_height,R0_height,R_height,PR0,PR,QT0,QT,ST0,ST] = Bxb_compare3(fbeats0, fbeats,pwaves0,pwaves,twaves0,twaves, acceptint, aOf, Of,swaves0,qwaves0,rwaves0)
% This function is similar to the function bxb.exe from Physionet's
% WFDB-App toolbox. It compares in a beat-by-beat basis if the detections
% match the reference. The algorithm is based on the entry by Joachim Behar
% on the Physionet / Computing in Cardiology Challenge 2013 and on ANSI/AAMI
% EC57 Norm 1998
%
% Input
% refqrs:        reference QRS detections
% testqrs:       detections to be tested against
% acceptint:     acceptance interval (left and right) in samples
%
%
% Output
% F1:            F1-measure (Joachim Behar - Computing in Cardiology 2013)
% ACC:           accuracy (by Karvounis 2007) - alternative to F1
% PPV:           positive predictive value
% SE:            sensitivity
% TP:            number of true positives
% FN:            number of false negatives
% FP:            number of false positives
%
% References
% [ANSI/AAMI Norm]  American National Standard ANSI/AAMI EC57:1998, Testing and Reporting Performance
% Results of Cardiac Rhythm and ST Segment Measurement Algorithms
%
% [WFDB-APP] Silva, I, Moody, G. "An Open-source Toolbox for Analysing and Processing PhysioNet Databases
% in MATLAB and Octave." Journal of Open Research Software 2(1):e27 [http://dx.doi.org/10.5334/jors.bi];
% 2014 (September 24).
%
% [Behar2014] Behar, J., Oster, J., & Clifford, G. D. (2014). Combining and Benchmarking Methods of Foetal
% ECG Extraction Without Maternal or Scalp Electrode Data. Physiological Measurement, 35(8), 1569??589.
%
%
% --
% fecgsyn toolbox, version 1.1, March 2016
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
%
% For more information visit: https://www.physionet.org/physiotools/ipmcode/fecgsyn/
%
% Referencing this work
%
%   Behar Joachim, Andreotti Fernando, Zaunseder Sebastian, Li Qiao, Oster Julien, Clifford Gari D.
%   An ECG simulator for generating maternal-foetal activity mixtures on abdominal ECG recordings.
%   Physiological Measurement.35 1537-1550. 2014.
%
% Last updated : 10-03-2016
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% == input test

if size(fbeats0,2) > size(fbeats0,1)
    fbeats0 = fbeats0';
end

if size(fbeats,2) > size(fbeats,1)
    fbeats = fbeats';
end

if size(pwaves0,2) > size(pwaves0,1)
    pwaves0 = pwaves0';
end

if size(pwaves,2) > size(pwaves,1)
    pwaves = pwaves';
end


if size(twaves0,2) > size(twaves0,1)
    twaves0 = twaves0';
end

if size(twaves,2) > size(twaves,1)
    twaves = twaves';
end

if size(rwaves0,2) > size(rwaves0,1)
    rwaves0 = rwaves0';
end

if size(swaves0,2) > size(swaves0,1)
    swaves0 = swaves0';
end
if size(qwaves0,2) > size(qwaves0,1)
    qwaves0 = qwaves0';
end


NB_REF = length(fbeats0);
NB_TEST = length(fbeats);

% == core function

[idxmatch,dist] = dsearchn(fbeats0,fbeats);     % closest ref for each point in test qrs
%amp_dist = aOf(refqrs(idxmatch)) - Of(testqrs);
fbeats0 = fbeats0(idxmatch(dist<acceptint));         % keep only the ones within a certain window
fbeats = fbeats(dist<acceptint);         % keep only the ones within a certain window

pw0_ind = [];
pw_ind = [];
tw0_ind = [];
tw_ind = [];
sw0_ind = [];
qw0_ind = [];

r_ind = [];
for i = 1:length(fbeats)
    if sum(max(find(pwaves0<fbeats0(i)))) == 0 || sum(max(find(pwaves<fbeats(i)))) == 0 ||sum(min(find(twaves0>fbeats0(i)))) == 0||sum(min(find(twaves>fbeats(i)))) == 0
        continue;
    end
    %[~,ind] = min(abs(rwaves0-fbeats(i)));
    r_ind = [r_ind;i];
    [~,ind] = min(abs(pwaves0-fbeats0(i)));
    pw0_ind =[pw0_ind;ind];
    [~,ind] = min(abs(pwaves-fbeats(i)));
    pw_ind =[pw_ind;ind];
    [~,ind] = min(abs(twaves0-fbeats0(i)));
    tw0_ind =[tw0_ind;ind];
    [~,ind] = min(abs(twaves-fbeats(i)));
    tw_ind =[tw_ind;ind];
    [~,ind] = min(abs(swaves0-fbeats0(i)));
    sw0_ind =[sw0_ind;ind];
    [~,ind] = min(abs(qwaves0-fbeats0(i)));
    qw0_ind =[qw0_ind;ind];
    
end

pw0 = pwaves0(pw0_ind);
tw0 = twaves0(tw0_ind);
pw = pwaves(pw_ind);
tw = twaves(tw_ind);
sw0 = swaves0(sw0_ind);
sw = sw0;
qw0 = qwaves0(qw0_ind);
qw = qw0;
fbeats0 = fbeats0(r_ind);
fbeats = fbeats(r_ind);
PR0 = fbeats0 - pw0;
PR = fbeats0 - pw;
QT0 = tw0-qw0;
QT = tw-qw;
ST0 = tw0-sw0;
ST = tw-sw;
isoe0 = round((tw0(1:(end-1))+pw0(2:end))/2);

%isoe = round((tw(1:(end-1))+pw(2:end))/2);

T0_height = aOf(tw0)-median(aOf(isoe0));
T_height = Of(tw)-median(aOf(isoe0));
P0_height = aOf(pw0)-median(aOf(isoe0));
P_height = Of(pw)-median(aOf(isoe0));
%T0_height = aOf(tw0);
%T_height = Of(tw);
%P0_height = aOf(pw0);
%P_height = Of(pw);
R0_height = abs(aOf(fbeats)-median(aOf(isoe0)));
R_height = abs(Of(fbeats)-median(aOf(isoe0)));
FTQRS0 = abs(T0_height./R0_height);
FTQRS = abs(T_height./R_height);
FTQRS_DR = sum(abs(FTQRS-FTQRS0)./abs(FTQRS0))/max(size(FTQRS0));


FPQRS0 = P0_height./R0_height;
FPQRS = P_height./R_height;
FPQRS_DR = sum(abs(FPQRS-FPQRS0)./abs(FPQRS0))/max(size(FPQRS0));



PR_DR = sum(abs((PR0-PR)./PR0))/max(size(PR0));
QT_DR = sum(abs((QT0-QT)./QT0))/max(size(QT0));
ST_DR = sum(abs((ST0-ST)./ST0))/max(size(ST0));
PR0 = PR0';
PR = PR';
QT0 = QT0';
QT = QT';
ST0 = ST0';
ST = ST';