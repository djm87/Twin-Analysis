%% Initialize Matlab
tic
clear all
close all

%% EBSD import convention
% Crystal symmetry
CS = {'notIndexed',crystalSymmetry('622', [2.95 2.95 4.6855],...
    'X||a', 'Y||b*', 'Z||c', 'mineral', 'Ti', 'color', 'light blue')};

% EBSD import 
dataDir='E:\Dropbox\Research\Source\Twin-Analysis\Data\HP-Ti'
% EBSDinterface='mat';
EBSDinterface='osc';
% EBSDinterface='ang';
EBSDfname = fullfile(dataDir,strcat('800C_5pct_TD_Cleaned.',EBSDinterface));
% EBSDfname = 'E:\Box Sync\PROJECT - Ti Compression\Clean_and_Rotate_EBSD\rotatedCleaned\RD_625C_5_rebsd.mat'

EBSDframe='convertSpatial2EulerReferenceFrame';

% EBSD frame
setMTEXpref('xAxisDirection','south');
setMTEXpref('zAxisDirection','outofPlane');


%% Cleanup Parameters
% CI threshold 
min_CI=0.08;

% Minimum grain size
min_grainSz = 1.5; %minimum indexed points per grain 

%% Edges to ignore 
edgeList2Remove=[]; % number of edge from label 

%% Grain Reconstruction Parameters
% grain boundary reconstruction
seg_angle = 5*degree; %misorientation of for determining fragments
seg_angle_grouped = 15*degree; % misorientation for grouping fragments (just family

%Misorientation tolerence for ebsd pixels being grouped
Mistol=5*degree; %<---- CHECK how this is used!

%% Twin boundary tolerances
meanMistol=8*degree; %needs to be slight larger (like Rod noted)
meanMistolRelaxed=10*degree; %used when including twin boundaries


%% Twin definitions
%Notes: All HCP are compound twins. Don't both with the types.
%Specify the twins
 twin={};
 tnum=1;
 twin{tnum}.CS=CS{2};
 twin{tnum}.name='T1 <10-1-1>(10-12)';
 twin{tnum}.k1=Miller(1,0,-1,2,twin{tnum}.CS,'hkl'); %What you specify here affects sign and schmid value
 twin{tnum}.eta1=Miller(-1,0,1,1,twin{tnum}.CS,'uvtw'); %What you specify here affects sign and schmid value
 twin{tnum}.actType=1; %i.e. 180 around K1, 180 around eta, compound K1 then eta1, compound eta1 then K1
 twin{tnum}.variantsToUse=1 

 tnum=2;
 twin{tnum}.CS=CS{2};
 twin{tnum}.name='T2 <-1-126>(11-21)';
 twin{tnum}.k1=Miller(-1,-1,2,1,twin{tnum}.CS,'hkl'); %What you specify here affects sign and schmid value
 twin{tnum}.eta1=Miller(1,1,-2,6,twin{tnum}.CS,'uvtw'); %What you specify here affects sign and schmid value
 twin{tnum}.actType=1; %i.e. 180 around K1, 180 around eta, compound K1 then eta1, compound eta1 then K1
 twin{tnum}.variantsToUse=1 
 
 tnum=3;
 twin{tnum}.CS=CS{2};
 twin{tnum}.name='C1 <11-2-3>(11-22)';
 twin{tnum}.k1=Miller(-1,-1,2,2,twin{tnum}.CS,'hkl'); %What you specify here affects sign and schmid value
 twin{tnum}.eta1=Miller(-1,-1,2,-3,twin{tnum}.CS,'uvtw'); %What you specify here affects sign and schmid value
 twin{tnum}.actType=1; %(pick 1-4):180 around K1, 180 around eta1, compound - K1 then eta1, compound - eta1 then K1
  twin{tnum}.variantsToUse=1 
  
 tnum=4;
 twin{tnum}.CS=CS{2};
 twin{tnum}.name='C2 <10-12>(10-1-1)';
 twin{tnum}.k1=Miller(-1,0,1,1,twin{tnum}.CS,'hkl'); %What you specify here affects sign and schmid value
 twin{tnum}.eta1=Miller(1,0,-1,2,twin{tnum}.CS,'uvtw'); %What you specify here affects sign and schmid value
 twin{tnum}.actType=1; %(pick 1-4):180 around K1, 180 around eta1, compound - K1 then eta1, compound - eta1 then K1
 twin{tnum}.variantsToUse=1 
 
 %Compute twin properties 
 twin=getTwinProperties(twin);
 
 %Visualize twins convention as described in Li eta al. 2015 
 %
%  tnum=3;
%  figure(19);plot(crystalShape.hex(CS{2}),'FaceAlpha',0.5)
%  hold on 
%  arrow3d(1.0*normalize(vector3d(twin{tnum}.k1)),'facecolor','black')
%  arrow3d(0.8*normalize(vector3d(twin{tnum}.eta1)),'facecolor','black')
%  arrow3d(normalize(vector3d(twin{tnum}.Rtw * twin{tnum}.k1)),'facecolor','green')
%  arrow3d(normalize(vector3d(twin{tnum}.Rtw * twin{tnum}.eta1)),'facecolor','green')
%  arrow3d(0.8*normalize(vector3d(CS{2}.aAxis)),'facecolor','red')
% %  arrow3d(0.8*normalize(vector3d(CS{2}.bAxis)),'facecolor','red')
%  arrow3d(0.8*normalize(vector3d(CS{2}.cAxis)),'facecolor','red')
%  hold off
%

%% Specimen stress state
sigma = stressTensor([0 0 0; 0 0 0; 0 0 -1]) %Sign of loading does more than just invert twin/parent flag