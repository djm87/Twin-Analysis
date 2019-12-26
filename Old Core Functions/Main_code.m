%% Initialize Matlab
tic
clear all
close all
%% Initialize run parameters 

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                       Set loadEBSD parameters
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% Crystal symmetry
CS = {'notIndexed',crystalSymmetry('622', [2.95 2.95 4.6855],...
    'X||a', 'Y||b*', 'Z||c', 'mineral', 'Ti', 'color', 'light blue')};

% EBSD import 
filename = 'Zr_LN_5pct_06'
twinAnalysisDatabase='E:\Dropbox\Research\Source\Twin-Analysis\Database\';
EBSDfname = ['E:\Dropbox\Research\Source\Twin-Analysis\Rods Data\Zr-LN-IP-5\LN-IP-05-6\' filename '.ang'];
% EBSDfname = 'E:\Dropbox\Research\Source\Twin-Analysis\Example Data\625C-RD-20pct cleaned-Dil cropped.ang'
EBSDinterface='ang';
% EBSDinterface='osc';
EBSDframe='convertSpatial2EulerReferenceFrame';

% EBSD frame
setMTEXpref('xAxisDirection','south');
setMTEXpref('zAxisDirection','outofPlane');
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                       Set Cleanup Parameters
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% CI threshold 
min_CI=0.08;

% Minimum grain size
min_grainSz = 1.5; %minimum indexed points per grain 
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                  Set Grain Reconstruction Parameters
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
% grain boundary reconstruction
seg_angle = 5*degree; %misorientation of for determining fragments
seg_angle_grouped = 15*degree; % misorientation for grouping fragments (just family


%Misorientation tolerence for ebsd pixels being grouped
Mistol=5*degree; %<---- CHECK how this is used!
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%                  Set Twin Boundary Parameters
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%Misorientation tolerence for mean orientations
meanMistol=8*degree; %needs to be slight larger (like Rod noted)
meanMistolRelaxed=15*degree; %used when including twin boundaries
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%
%Notes: All HCP are compound twins. Don't both with the types.
%Specify the twins
 twin={};
 tnum=1;
 twin{tnum}.CS=CS{2};
 twin{tnum}.name='T1 <10-1-1>(10-12)';
 twin{tnum}.k1=Miller(1,0,-1,2,twin{tnum}.CS,'hkl'); %What you specify here affects sign and schmid value
 twin{tnum}.eta1=Miller(-1,0,1,1,twin{tnum}.CS,'uvtw'); %What you specify here affects sign and schmid value
 twin{tnum}.actType=1; %i.e. 180 around K1, 180 around eta, compound K1 then eta1, compound eta1 then K1
 twin{tnum}.variantsToUse=1; 

 tnum=2;
 twin{tnum}.CS=CS{2};
 twin{tnum}.name='T2 <-1-126>(11-21)';
 twin{tnum}.k1=Miller(-1,-1,2,1,twin{tnum}.CS,'hkl'); %What you specify here affects sign and schmid value
 twin{tnum}.eta1=Miller(1,1,-2,6,twin{tnum}.CS,'uvtw'); %What you specify here affects sign and schmid value
 twin{tnum}.actType=1; %i.e. 180 around K1, 180 around eta, compound K1 then eta1, compound eta1 then K1
 twin{tnum}.variantsToUse=1;
 
 tnum=3;
 twin{tnum}.CS=CS{2};
 twin{tnum}.name='C1 <11-2-3>(11-22)';
 twin{tnum}.k1=Miller(-1,-1,2,2,twin{tnum}.CS,'hkl'); %What you specify here affects sign and schmid value
 twin{tnum}.eta1=Miller(-1,-1,2,-3,twin{tnum}.CS,'uvtw'); %What you specify here affects sign and schmid value
 twin{tnum}.actType=1; %(pick 1-4):180 around K1, 180 around eta1, compound - K1 then eta1, compound - eta1 then K1
  twin{tnum}.variantsToUse=1;
 
%  tnum=4;
%  twin{tnum}.CS=CS{2};
%  twin{tnum}.name='C2 <10-12>(10-1-1)';
%  twin{tnum}.k1=Miller(-1,0,1,1,twin{tnum}.CS,'hkl'); %What you specify here affects sign and schmid value
%  twin{tnum}.eta1=Miller(1,0,-1,2,twin{tnum}.CS,'uvtw'); %What you specify here affects sign and schmid value
%  twin{tnum}.actType=1; %(pick 1-4):180 around K1, 180 around eta1, compound - K1 then eta1, compound - eta1 then K1
%   twin{tnum}.variantsToUse=1 
 %Compute twin properties 
 twin=getTwinProperties(twin);
 
 %Visualize twins convention as described in Li eta al. 2015 
 %
%  tnum=3;
%  figure(19);plot(crystalShape.hex(CS{2}),'FaceAlpha',0.5)
%  hold on 
%  arrow3d(1.2*normalize(vector3d(twin{tnum}.k1)),'facecolor','black')
%  arrow3d(0.8*normalize(vector3d(twin{tnum}.eta1)),'facecolor','black')
%  arrow3d(normalize(vector3d(twin{tnum}.Rtw * twin{tnum}.k1)),'facecolor','green')
%  arrow3d(normalize(vector3d(twin{tnum}.Rtw * twin{tnum}.eta1)),'facecolor','green')
%  arrow3d(0.8*normalize(vector3d(CS{2}.aAxis)),'facecolor','red')
% %  arrow3d(0.8*normalize(vector3d(CS{2}.bAxis)),'facecolor','red')
%  arrow3d(0.8*normalize(vector3d(CS{2}.cAxis)),'facecolor','red')
%  hold off
 %
 %Specify specimen stress state
sigma = stressTensor([0 0 0; 0 -1 0; 0 0 0]) %Sign of loading does more than just invert twin/parent flag

%% Import the Data
% create an EBSD variable containing the data
ebsd = loadEBSD(EBSDfname,CS,'interface',EBSDinterface,EBSDframe);
% ebsd=ebsd.gridify;
%% Clean Data
%color scheme of ipf maps
oM = ipdfHSVOrientationMapping(ebsd('indexed'));
% figure;plot(oM)
color = oM.orientation2color(ebsd('indexed').orientations);

%Raw data 
% figure; plot(ebsd,ebsd.orientations)

%Cleaned data
ebsd(ebsd.prop.iq==0).phaseId=1;
ebsd=ebsd(ebsd.prop.ci>min_CI);
figure; plot(ebsd,ebsd.orientations)

%% Reconstruct and clean grains
%Calculate all grains
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),...
    'angle',seg_angle,'boundary','tight');

% figure;plot(grains('Mg'),grains('Mg').meanOrientation)

%Remove grains that are too small
grainsSelected = grains(grains.grainSize > min_grainSz);

%reindex ebsd 
ebsd = ebsd(grainsSelected);

% and perform grain reconstruction with the reduced EBSD data set
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'angle',seg_angle);

% figure;plot(grains('Mg'),grains('Mg').meanOrientation)

%% Cluster grains based on fragment misorientation and compute grain fragment properties
% Use a graph to handle and visualize the data. The motivation for using a
% graph is data management and visualization.

%Takes about 1 minute with 1500 grains
doplot=true;

G = InitializeGraph(ebsd,grains,twin,Mistol,meanMistol,meanMistolRelaxed,doplot)

%% Create user list to delete incorrect twin relationships
%Remove edges that that don't look right for a twin.
edgeList2Remove=[]; % number of edge from label 

doplot=true;
dolabel=true;
G_clust = ClusterGrainsTwins(G,grains,Mistol,meanMistolRelaxed,...
    twin,edgeList2Remove,doplot,dolabel)

%Call first time 
w=[1,1,1];

G_tmp = BuildTwinRelationships(G_clust,grains,w,twin,...
    sigma,seg_angle_grouped,Mistol,meanMistolRelaxed);

%% Extract parent stats for schmid and compute the angle
G_Complete = GetSchmidVariants(G_tmp,twin,sigma);

%% Reconstruct EBSD
[~,G_Complete] = fragmentReconstruction(G_Complete,ebsd,grains);
figure; plot(grains,G_Complete.Nodes.MergeTwin)

%% Get current count statistics 
[stats.twinCount,G_Complete] = CountTwins(G_Complete)
% figure; 
% hist(stats.twinCount)
% figure; 
% plot(grains,G_Complete.Nodes.twinCount)

%% Get twin thickness 
[G_Complete] = TwinThickness(G_Complete,grains)

selection=G_Complete.Nodes.twinThickness<1;
figure; [~,mP] = plot(grains(selection),G_Complete.Nodes.twinThickness(selection))
mP.micronBar.length=2;
%% Compute the volume fractions 
[stats.twinVF,stats.totalArea,stats.totalTwinVF] = ...
    getTwinFractions(G_Complete,grains,twin)

%% Transfer G to grains
[grains] = transferGtoGrains(G_Complete,grains)
%%
figure; plot(grains,grains.prop.twinCount)
%% Save files for processing later 
save([twinAnalysisDatabase filename '_twin_database1.mat'],'grains','stats','twin')
% %% Make some pretty plots
% [mergedGrains] = MergeByEdge(G_Complete,ones(length(G_Complete.Edges.pairs),1,'logical'),grains)
% 
% figure; 
% plot(grains,G_Complete.Nodes.Type,'Micronbar','off')
% hold on 
% plot(mergedGrains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-')
% hold off
% mtexColorbar
% figure; 
% plot(grains,G_Complete.Nodes.Generation,'Micronbar','off')
% hold on 
% plot(mergedGrains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-')
% hold off
% mtexColorbar
% figure; 
% plot(grains,grains.meanOrientation,'Micronbar','off')
% hold on 
% plot(mergedGrains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-')
% hold off
% 
% figure;
% hist(G_Complete.Edges.schmidActive) %Value of active schmid
% figure;
% hist(G_Complete.Edges.schmidActiveN) %variant of rotation of the axis variant that is active
% figure;
% hist(G_Complete.Edges.schmidActiveRank) %variant order in decreasing SF
% figure; 
% plot(grains,G_Complete.Nodes.meanOrientation ,'Micronbar','off')
% hold on 
%  p=plot(G_Complete,'XData',G_Complete.Nodes.centroids(:,1),...
%             'YData',G_Complete.Nodes.centroids(:,2),'displayName','graph');
% p.Marker='s';p.NodeColor='k';p.MarkerSize=3;p.EdgeColor='k';
% labeledge(p,G_Complete.Edges.pairs(:,1),...
%     G_Complete.Edges.pairs(:,2),G_Complete.Edges.GlobalID);
% labelnode(p,G_Complete.Nodes.Id,G_Complete.Nodes.Id);
% hold off
% 
% figure; 
% plot(grains,G_clust.Nodes.meanOrientation ,'Micronbar','off')
% hold on 
%  p=plot(G_clust,'XData',G_clust.Nodes.centroids(:,1),...
%             'YData',G_clust.Nodes.centroids(:,2),'displayName','graph');
% p.Marker='s';p.NodeColor='k';p.MarkerSize=3;p.EdgeColor='k';
% labeledge(p,G_clust.Edges.pairs(:,1),...
%     G_clust.Edges.pairs(:,2),1:length(G_clust.Edges.pairs));
% labelnode(p,G_clust.Nodes.Id,G_clust.Nodes.Id);
% hold off
% figure;plot(grains2,grains2.meanOrientation)
% figure;plot(grains,grains.meanOrientation)

% figure; plot(grains(G_Complete.Nodes.MergeTwin),grains(G_Complete.Nodes.MergeTwin).meanOrientation)