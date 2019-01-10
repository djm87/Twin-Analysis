%% Initialize Matlab
clear all
close all
%% Initialize run parameters 
%takes about 2 minutes for 1500 grains
tic
% plotting conventions
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outofPlane');

% ebsd import 
EBSDfname = 'Area3.ang';
EBSDinterface='ang';
EBSDframe='convertEuler2SpatialReferenceFrame';

% crystal symmetry
% CS = {crystalSymmetry('622', [3.2 3.2 5.2], 'X||a', 'Y||b*', 'Z||c', 'mineral', 'Mg', 'color', 'light blue')};
CS = {'notIndexed',...
  crystalSymmetry('622', [3.2 3.2 5.2], 'X||a', 'Y||b*', 'Z||c', 'mineral', 'Mg', 'color', 'light blue')};
% CI threshold 
min_CI=0.08;

% grain boundary reconstruction
seg_angle = 5*degree; %misorientation of for determining fragments
seg_angle_grouped = 15*degree; % misorientation for same fragments
min_grainSz = 1.5; %minimum indexed points per grain 

%Misorientation tolerence for ebsd pixels 
Mistol=5*degree;
%Misorientation tolerence for mean orientations
meanMistol=8*degree; %needs to be slight larger (like Rod noted)
meanMistolRelaxed=15*degree; %used when including twin boundaries

%Specify the twins
twins = {orientation('axis',Miller(1,1,-2,0,CS{2},'uvw'),...
        'angle',86.3471*degree,CS{2},CS{2}),... %Extension twin <10-11> largest amount
   orientation('axis',Miller(1,0,-1,0,CS{2},'uvw'),...
        'angle',34.7*degree,CS{2},CS{2}),... %Extension twin <1126> 
   orientation('axis',Miller(1,1,-2,0,CS{2},'uvw'),...
        'angle',56.0*degree,CS{2},CS{2}),... %Compression
   orientation('axis',Miller(1,1,-2,0,CS{2},'uvw'),...
        'angle',37.5*degree,CS{2},CS{2}),... %Double A
   orientation('axis',Miller(1,1,-2,0,CS{2},'uvw'),...
        'angle',30.1*degree,CS{2},CS{2})}; %Double B
    
twinNames={'Tension 1';'Tension 2';'Compression 1';'Double 1';'Double 2';'Double 3';'Double 4';'Double 5';'Double 6'}
% Specify CRSS, slip systems and stress state for schmid computations
CRSS=[111,133,125,125,111,111,111,111];
% CRSS=[25,80,210];
% sS={slipSystem(Miller(1,1,-2,0,CS{2},'uvtw'), Miller(0,0,0,1,CS{2},'hkl'),CRSS(1)),...
%     slipSystem(Miller(2,-1,-1,0,CS{2},'uvtw'), Miller(0,1,-1,0,CS{2},'hkl'),CRSS(2)),...
%     slipSystem(Miller(2,-1,-1,3,CS{2},'uvtw'), Miller(-1,1,0,1,CS{2},'hkl'),CRSS(3))};
sS={slipSystem(Miller(1,0,-1,-1,CS{2},'uvtw'), Miller(1,0,-1,2,CS{2},'hkl'),CRSS(1)),...
    slipSystem(Miller(-1,-1,2,6,CS{2},'uvtw'), Miller(1,1,-2,1,CS{2},'hkl'),CRSS(2)),...
    slipSystem(Miller(1,0,-1,-2,CS{2},'uvtw'), Miller(1,0,-1,1,CS{2},'hkl'),CRSS(3)),...
    slipSystem(Miller(1,0,-1,-1,CS{2},'uvtw'), Miller(1,0,-1,2,CS{2},'hkl'),CRSS(1)),...
    slipSystem(Miller(1,0,-1,-1,CS{2},'uvtw'), Miller(1,0,-1,2,CS{2},'hkl'),CRSS(1))};

%apply symmetries for full slip systems
for i=1:length(sS) 
    sS{i}=sS{i}.symmetrise('antipodal');
end
sS{3}
% sigma = tensor([0 0 0; 0 0 0; 0 0 -1],'name','stress') %Sign of loading does more than just invert twin/parent flag
sigma = stressTensor.uniaxial(-vector3d.Z)
%% Import the Data
% create an EBSD variable containing the data
ebsd = loadEBSD(EBSDfname,CS,'interface',EBSDinterface,EBSDframe);
%% Clean Data
%color scheme of ipf maps
oM = ipdfHSVOrientationMapping(ebsd('indexed'));
figure;plot(oM)
color = oM.orientation2color(ebsd('indexed').orientations);

%Raw data 
figure; plot(ebsd,ebsd.orientations)

%Cleaned data
ebsd(ebsd.prop.iq==0).phaseId=1;
ebsd=ebsd(ebsd.prop.ci>min_CI);
figure; plot(ebsd,ebsd.orientations)
%% Reconstruct and clean grains
%Calculate all grains
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),...
    'angle',seg_angle,'boundary','tight');

figure;plot(grains('Mg'),grains('Mg').meanOrientation)

%Remove grains that are too small
grainsSelected = grains(grains.grainSize > min_grainSz);

%reindex ebsd 
ebsd = ebsd(grainsSelected);

% and perform grain reconstruction with the reduced EBSD data set
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'angle',seg_angle);

figure;plot(grains('Mg'),grains('Mg').meanOrientation)
%% Cluster grains based on fragment misorientation and compute grain fragment properties
% Use a graph to handle and visualize the data. The motivation for using a
% graph is data management and visualization.

%Takes about 1 minute with 1500 grains
doplot=true;

G = InitializeGraph(ebsd,grains,twins,Mistol,meanMistol,meanMistolRelaxed,doplot)

%% Cluster grains based on twin boundaries alone

doplot=true;
dolabel=true;

G_clust_twin = ClusterGrainsTwins(G,grains,Mistol,meanMistolRelaxed,...
    twins,[],doplot,dolabel)

%% Create user list to delete incorrect twin relationships
%Remove edges that that don't look right for a twin.
edgeList2Remove=[785;764;1461;568;571;572;1727;2004;1770;1552;]; % number of edge from label 

doplot=true;
dolabel=true;
G_clust = ClusterGrainsTwins(G,grains,Mistol,meanMistolRelaxed,...
    twins,edgeList2Remove,doplot,dolabel)
%% To Do: Add twin edges
% TO DO:This is straightforward. However, getting twin type right requires some
% sort of interactive misorientation analysis which is not obvious how to
% do in MTEX. OIM is better for this. 
%% Compute Schmid info for twin/parents in clustered grains
%This computes Schmid factor for twin/parent identification

G_clust = GetSchmidRelative(G_clust,twins,CRSS,sS,sigma);

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
w=[0,1,3]
G_Complete = FamilyVotes(G_clust,w);


%% Make the family tree
G_Complete = MakeFamilyTree(G_Complete);

% Compare results with boundaries 
% Mistol=5*degree
[twinBoundary] = GetTwinBoundaries(G_Complete,grains,twins,twinNames,8*degree);
% toc/60
% Need to rethink how grains are selected 

%% Compute the volume fractions 
Areas=G_Complete.Nodes.Area;
T1=sum(Areas(G_Complete.Nodes.Type==1))/ sum(area(grains)) * 100
T2=sum(Areas(G_Complete.Nodes.Type==2))/ sum(area(grains)) * 100
C=sum(Areas(G_Complete.Nodes.Type==3))/ sum(area(grains)) * 100
D1=sum(Areas(G_Complete.Nodes.Type==4))/ sum(area(grains)) * 100
D2=sum(Areas(G_Complete.Nodes.Type==5))/ sum(area(grains)) * 100
[T1;T2;C;D1;D2]
sum(area(grains))



