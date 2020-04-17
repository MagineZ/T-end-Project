This directory contains Matlab files for R peak detection and T-wave end location,
as well as demo files.

The most important files are rpeak.m for R peak detection and twaveend.m for
T-wave end location.

To run the demos, the PhysioNet QT database converted into text format
(located in the qtdata directory included in the same package) must be
available. 

IMPORTANT: Before running the demo files qtdemos.m and otherdemos.m in Matlab,
the file qtdatapath.m must be edited to indicate the absolute path of the 
qtdata directory (do not omit the file separator at the end of the path).

To reproduce the results reported in PI-1744.pdf, run the Matlab files blpsample1.m,
blpsample2.m, blprecord1.m and blprecord2.m.

Author: Qinghua Zhang
Copyright 2005 INRIA

