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
CS = {'notIndexed',crystalSymmetry('622', [3.224 3.224 5.212],...
    'X||a', 'Y||b*', 'Z||c', 'mineral', 'WE43', 'color', 'light blue')};

% EBSD import 
EBSDfname = 'E:\Dropbox\Research\WE43\TwinAnalysis\Analyzed\03-WE43-T6- Actapaper-high\TD\0p15.osc';
EBSDinterface='osc';
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
seg_angle_grouped = 15*degree; % misorientation for grouping fragments


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
 twin{tnum}.name='C1 (10-12)<10-1-1>';
 twin{tnum}.k1=Miller(-1,0,1,1,twin{tnum}.CS,'hkl'); %What you specify here affects sign and schmid value
 twin{tnum}.eta1=Miller(1,0,-1,2,twin{tnum}.CS,'uvtw'); %What you specify here affects sign and schmid value
 twin{tnum}.actType=1; %i.e. 180 around K1, 180 around eta, compound K1 then eta1, compound eta1 then K1
 twin{tnum}.variantsToUse=1;
 
 tnum=4;
 twin{tnum}.CS=CS{2};
 twin{tnum}.name='DT (10-12)-(10-1-1)';
 twin{tnum}.k1(2)=Miller(1,0,-1,2,twin{tnum}.CS,'hkl'); %What you specify here affects sign and schmid value
 twin{tnum}.k1(1)=Miller(-1,0,1,1,twin{tnum}.CS,'hkl'); %What you specify here affects sign and schmid value
 twin{tnum}.eta1(2)=Miller(-1,0,1,1,twin{tnum}.CS,'uvtw'); %What you specify here affects sign and schmid value
 twin{tnum}.eta1(1)=Miller(1,0,-1,2,twin{tnum}.CS,'uvtw'); %What you specify here affects sign and schmid value
 twin{tnum}.actType=1; %i.e. 180 around K1, 180 around eta   
 twin{tnum}.variantsToUse=[1,2];

 %Compute twin properties 
 twin=getTwinProperties(twin);
  
 %Visualize twins convention as described in Li eta al. 2015 
 tnum=1;
 figure;plot(crystalShape.hex(CS{2}),'FaceAlpha',0.5)
 hold on 
 arrow3d(normalize(vector3d(twin{tnum}.k1)),'facecolor','black')
 arrow3d(normalize(vector3d(twin{tnum}.eta1)),'facecolor','black')
 arrow3d(normalize(vector3d(twin{tnum}.Rtw * twin{tnum}.k1)),'facecolor','green')
 arrow3d(normalize(vector3d(twin{tnum}.Rtw * twin{tnum}.eta1)),'facecolor','green')
 arrow3d(0.8*normalize(vector3d(CS{2}.aAxis)),'facecolor','red')
%  arrow3d(0.8*normalize(vector3d(CS{2}.bAxis)),'facecolor','red')
 arrow3d(0.8*normalize(vector3d(CS{2}.cAxis)),'facecolor','red')
 hold off

 %Specify specimen stress state
 sigma = tensor([0 0 0; 0 0 0; 0 0 -1],'name','stress') %Sign of loading does more than just invert twin/parent flag
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
ebsd(ebsd.prop.imagequality==0).phaseId=1;
ebsd=ebsd(ebsd.prop.confidenceindex>min_CI);
% figure; plot(ebsd,ebsd.orientations)

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

%% Cluster grains based on twin boundaries alone

doplot=true;
dolabel=true;

G_clust_twin = ClusterGrainsTwins(G,grains,Mistol,meanMistolRelaxed,...
    twin,[],doplot,dolabel)

%% Create user list to delete incorrect twin relationships
%Remove edges that that don't look right for a twin.
edgeList2Remove=[577,543,5643,7827,7814,2269,1288,]; % number of edge from label 

doplot=true;
dolabel=true;
G_clust = ClusterGrainsTwins(G,grains,Mistol,meanMistolRelaxed,...
    twin,edgeList2Remove,doplot,dolabel)

%% Compute Schmid info for twin/parents in clustered grains
%This computes Schmid factor for twin/parent identification
G_clust = GetSchmidRelative(G_clust,twin,sigma);

tmp=G_clust.Nodes.EffSF(:,1);
notZero=tmp<1;
  figure; 
    plot(grains(G_clust.Nodes.Id(notZero)),...
        tmp(notZero),'Micronbar','off');
    hold on 
    p=plot(G_clust,'XData',G_clust.Nodes.centroids(:,1),...
        'YData',G_clust.Nodes.centroids(:,2));
    hold off

%% Identify Families in grain clusters 
%Remove nodes that aren't connected by edges, alternatively use condition
%on minimum occurance from conncomp
doplot=true;
dolabel=false;
G_clust = AssignFamilyIDs(G_clust,grains,seg_angle_grouped,doplot,dolabel);
%% Implement voting scheme for the respective families
% For each family compute the vote between connected families 
%Edges currently connect clustered grains. We can compute the vote for each
%edge and then average for the family. This ensures that only neighbors are
%ever voted on. So that votes are really between familys average values of
%schmid are computed and stored for each fragment in a family. The area and
%relative boundary are taken in the addative sense rather than average. 

%Edges.Gb contains boundaries between pairs 
%Nodes.Gb contains boundaries for each fragment 
%To compute the relative boundary between families we need to sum boundary
%of a particular type and since each pair is its own type this should be
%relatively simple. 
w=[1,1,1]
G_Complete_unclean = FamilyVotes(G_clust,w);


%% Cleanup the family tree
runCleanup=true
maxIter=5;
CleanupIter=0;
%I think this shouldn't take more than 2 itterations. There can be some
%cases that aren't solved until grains are regroupd (i.e you have two
%seperate grains that need still need to be split and a circular relationship exists
%because the same families are in both grains. It's a thing.. go figure!)

%  G_Complete_unclean=ClusterGrainsTwins(G_Complete_unclean,grains,Mistol,meanMistolRelaxed,...
%     twin,[438],false,false);

while runCleanup && CleanupIter<maxIter
    [G_Complete_unclean,runCleanup] = CleanFamilyTree(G_Complete_unclean,grains);
    
    G_Complete_unclean=ClusterGrainsTwins(G_Complete_unclean,grains,Mistol,meanMistolRelaxed,...
    twin,[],false,false);

    CleanupIter=CleanupIter+1;
end
if CleanupIter > maxIter
    disp('More CleanupIter than expected, likely need to solve something manually or improve the code')
    pause;
end
%% Reconstruct grains based on cleanup (some new grains may have been created
edgeList2Remove=[]; % number of edge from label 
% 92,86,216,361,299,392,224,235
doplot=true;
dolabel=true;
G_Complete = ClusterGrainsTwins(G_Complete_unclean,grains,Mistol,meanMistolRelaxed,...
    twin,edgeList2Remove,doplot,dolabel)
%% Redo Families in grain clusters 
%Remove nodes that aren't connected by edges, alternatively use condition
%on minimum occurance from conncomp
doplot=true;
dolabel=false;
G_Complete = AssignFamilyIDs(G_Complete,grains,seg_angle_grouped,doplot,dolabel);

%% Create the family tree 
%uses the recursive relationship as described in Pradalier et al. 2018
[G_Complete] = CreateFamilyTree(G_Complete,grains);
figure; 
plot(grains,G_Complete.Nodes.Type,'Micronbar','off')
%%
% Compare results with boundaries 
% Mistol=5*degree
[twinBoundary] = GetTwinBoundaries(G_Complete,grains,twin,8*degree);
% toc/60
% Need to rethink how grains are selected 

%% Compute the volume fractions 
Areas=G_Complete.Nodes.Area;
T1=sum(Areas(G_Complete.Nodes.Type==1))/ sum(area(grains)) * 100
T2=sum(Areas(G_Complete.Nodes.Type==2))/ sum(area(grains)) * 100
% C=sum(Areas(G_Complete.Nodes.Type==3))/ sum(area(grains)) * 100

sum([T1 T2])
toc


