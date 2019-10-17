%Use: in matlab or octave navigate to the git repository and execute this script. 
%This installs the git repository to the startup paths for Matlab or Octave. 
%All scripts can then be executed from any directory in you system, allowing you to
%avoid having to copy the scripts around for different refinement directories.

P=genpath(fullfile(pwd,'Core Functions'))
addpath(P)
savepath

run(fullfile(pwd,'MTEX','install_mtex'))