%% Initialize Matlab

clear all
close all

%% Set the parallel environment
opt.ncpu=4;
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolobj=parpool(opt.ncpu)
    poolobj.IdleTimeout=220;
end
%% EBSD import convention
% Crystal symmetry
opt.CS  = {'notIndexed',crystalSymmetry('622', [2.95 2.95 4.6855],...
    'X||a', 'Y||b*', 'Z||c', 'mineral', 'Ti', 'color', 'light blue')};
% EBSD import 
opt.EBSD.dataDir='E:\Dropbox\Research\Source\Twin-Analysis\Data\HP-Ti';
opt.EBSD.interface='osc';% 'osc','mat';'ang'
opt.EBSD.fname = fullfile(opt.EBSD.dataDir,strcat('800C_20pct_RD_Cleaned.',opt.EBSD.interface));
opt.EBSD.frame='convertSpatial2EulerReferenceFrame';

% EBSD plot preferences
setMTEXpref('xAxisDirection','south');
setMTEXpref('zAxisDirection','outofPlane');
setMTEXpref('generatingHelpMode',true);
% setMTEXpref('defaultColorMap',hsv)
setMTEXpref('showMicronBar','off')
%% Plot options 
opt.plot.do=true;
opt.plot.ClusterOnly=true;
opt.plot.labelEdges=false;
opt.plot.legendOn=false;
%% Grain Reconstruction Parameters

%misorientation of for determining fragments
opt.grain_recon.seg_angle = 4*degree; 

% misorientation for grouping fragments (just used for FamilyID)
opt.grain_recon.seg_angle_grouped = 15*degree; 

% CI threshold 
opt.grain_recon.min_CI=0.08;

% Minimum grain size
opt.grain_recon.min_grainSz = 1.5; %minimum indexed points per grain 
%% Twin definitions
%Notes: HCP twins are compound twins. actType doesn't matter
 twin={};
 tnum=1;
 twin{tnum}.CS=opt.CS{2};
 twin{tnum}.name='T1 <10-1-1>(10-12)';
 twin{tnum}.k1=Miller(1,0,-1,2,twin{tnum}.CS,'hkl'); %What you specify here affects sign and schmid value
 twin{tnum}.eta1=Miller(-1,0,1,1,twin{tnum}.CS,'uvtw'); %What you specify here affects sign and schmid value
 twin{tnum}.actType=1; %i.e. 180 around K1, 180 around eta, compound K1 then eta1, compound eta1 then K1
 twin{tnum}.variantsToUse=1; 
 twin{tnum}.tol.misGb=6*degree; %Used to produce the merged grains starting point for segmentation
 twin{tnum}.tol.misGbRlx=8*degree; %Used to produce the merged grains starting point for segmentation
 twin{tnum}.tol.GbMinLen=5; %if the gB with a twin relationship is less than GbMinLen then mean relationship is used instead
 twin{tnum}.tol.misMean=5*degree; %determine twins when gB length is less than GbMinLen
 twin{tnum}.tol.misMeanRlx=10*degree; %determines the type for all edge relations
 twin{tnum}.voteWeights=[1,0,0]; %Weights for schmid, boundary, and area  based voting

 tnum=2;
 twin{tnum}.CS=opt.CS{2};
 twin{tnum}.name='T2 <-1-126>(11-21)';
 twin{tnum}.k1=Miller(-1,-1,2,1,twin{tnum}.CS,'hkl'); %What you specify here affects sign and schmid value
 twin{tnum}.eta1=Miller(1,1,-2,6,twin{tnum}.CS,'uvtw'); %What you specify here affects sign and schmid value
 twin{tnum}.actType=1; %i.e. 180 around K1, 180 around eta, compound K1 then eta1, compound eta1 then K1
 twin{tnum}.variantsToUse=1; 
 twin{tnum}.tol.misGb=2*degree; %Used to produce the merged grains starting point for segmentation
 twin{tnum}.tol.misGbRlx=4*degree; %Used to produce the merged grains starting point for segmentation
 twin{tnum}.tol.misMean=5*degree; %determine twins when gB length is less than GbMinLen
 twin{tnum}.tol.misMeanRlx=8*degree; %determines the type for all edge relations
 twin{tnum}.voteWeights=[2,1,1]; %Weights for schmid, boundary, and area  based voting

 tnum=3;
 twin{tnum}.CS=opt.CS{2};
 twin{tnum}.name='C1 <11-2-3>(11-22)';
 twin{tnum}.k1=Miller(-1,-1,2,2,twin{tnum}.CS,'hkl'); %What you specify here affects sign and schmid value
 twin{tnum}.eta1=Miller(-1,-1,2,-3,twin{tnum}.CS,'uvtw'); %What you specify here affects sign and schmid value
 twin{tnum}.actType=1; %(pick 1-4):180 around K1, 180 around eta1, compound - K1 then eta1, compound - eta1 then K1
  twin{tnum}.variantsToUse=1; 
 twin{tnum}.tol.misGb=6*degree; %Used to produce the merged grains starting point for segmentation
 twin{tnum}.tol.misGbRlx=8*degree; %Used to produce the merged grains starting point for segmentation
 twin{tnum}.tol.misMean=5*degree; %determine twins when gB length is less than GbMinLen
 twin{tnum}.tol.misMeanRlx=10*degree; %determines the type for all edge relations
 twin{tnum}.voteWeights=[5,1,1]; %Weights for schmid, boundary, and area  based voting
 
 tnum=4;
 twin{tnum}.CS=opt.CS{2};
 twin{tnum}.name='C2 <10-12>(10-1-1)';
 twin{tnum}.k1=Miller(-1,0,1,1,twin{tnum}.CS,'hkl'); %What you specify here affects sign and schmid value
 twin{tnum}.eta1=Miller(1,0,-1,2,twin{tnum}.CS,'uvtw'); %What you specify here affects sign and schmid value
 twin{tnum}.actType=1; %(pick 1-4):180 around K1, 180 around eta1, compound - K1 then eta1, compound - eta1 then K1
 twin{tnum}.variantsToUse=1; 
 twin{tnum}.tol.misGb=2*degree; %Used to produce the merged grains starting point for segmentation
 twin{tnum}.tol.misGbRlx=4*degree; %Used to produce the merged grains starting point for segmentation
 twin{tnum}.tol.misMean=5*degree; %determines twins when gB length is less than GbMinLen
 twin{tnum}.tol.misMeanRlx=8*degree; %determines the type for all edge relations
 twin{tnum}.voteWeights=[2,1,1]; %Weights for schmid, boundary, and area  based voting
 
 %Compute twin properties 
 twin=getTwinProperties(twin);

 %store in opt and set the unknown twin type
 opt.twin=twin;
 opt.nTwin=length(twin)-1;
 opt.twinUnknown=length(twin);
 opt.GbMinLen=5;
 opt.useMeanSmallGb=true;
 opt.mergeInclusion=true; %in a grain frag
 opt.mergeInclusionCluster=true; %in a merged grain 
 opt.mergeTripplePoints=false;
%% Specimen stress state
 opt.sigma=stressTensor([0 0 0; 0 0 0; 0 0 -1]) %Sign of loading does more than just invert twin/parent flag

%% Cleanup local variables 
clear 'twin' 'tnum'
opt