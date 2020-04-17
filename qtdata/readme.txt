This directory contains the PhysioNet QT database converted into text format. 
The original data in binary format is available at http://www.physionet.org.

The original binary data have been converted into text format in order to help users who 
have difficulties to install the tools necessary for the use of the original binary data. 
Following the rules stated on the PhysioNet web site, the Matlab program qtdbconvert.m 
is included in this package for reproducing the data transformation.

In order to reproduce these data files, the commands rdann and rdsamp from
the PhysioNet WFDB package must be installed. It seems that WFDB works only 
on Linux or Unix computers. Following the following steps to reproduce the
data files:
- make sure that the commands rdann and rdsamp work correctly;
- download .dat, .hea, .q1c, .q2c files from the PhysioNet QT databas;
- place the downloaded files and the Matlab file qtdbconvert.m (included
  in this package) in the same directory;
- In Matlab, go to the directory containing the above files, and run qtdbconvert.m.

For each data file of the original PhysioNet QT database, only the portion
of the signal annotated by one or two cardilogists is converted to text format.




