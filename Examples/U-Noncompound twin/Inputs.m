%% Initialize Matlab
clear all
close all

%% Set the parallel environment
% opt.ncpu=4;
% poolobj = gcp('nocreate'); % If no pool, do not create new one.
% if isempty(poolobj)
%     poolobj=parpool(opt.ncpu)
%     poolobj.IdleTimeout=30;
% end
%% EBSD import convention
%Load the EBSD and reconstructed grain map
load('ebsd_recon.mat');
load('grains_recon.mat');

% Crystal symmetry
opt.CS=ebsd.CSList';

% EBSD plot preferences
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outofPlane');
setMTEXpref('generatingHelpMode',true);
setMTEXpref('showMicronBar','off')

%% Twin definitions and boundary merging parameters
%Notes: HCP twins are compound twins. actType doesn't matter
 twin={};
 tnum=1;
 twin{tnum}.CS=opt.CS{2};
 twin{tnum}.name='{130}';
 twin{tnum}.k1=Miller(1,3,0,twin{tnum}.CS,'hkl'); %What you specify here affects sign and schmid value
 twin{tnum}.eta1=Miller(3,-1,0,twin{tnum}.CS,'uvw'); %What you specify here affects sign and schmid value
 twin{tnum}.actType=1; %i.e. 180 around K1, 180 around eta, compound K1 then eta1, compound eta1 then K1
 twin{tnum}.variantsToUse=1; 
 twin{tnum}.tol.misGb=5*degree; %Used to produce the merged grains starting point for segmentation
 twin{tnum}.tol.misGbRlx=10*degree; %Used to 
 twin{tnum}.tol.misMean=10*degree; %determine twins when gB length is less than GbMinLen
 rot=twin{tnum}.CS.rotation
 twin{tnum}.variantRotations=rot(:,1);
 twin{tnum}.voteWeights=[2,1,1,1]; %Weights for schmid, boundary, area, and initial ODF
 
 tnum=2;
 twin{tnum}.CS=opt.CS{2};
 twin{tnum}.name='{172}';
 twin{tnum}.k1=Miller(1,7,2,twin{tnum}.CS,'hkl'); %What you specify here affects sign and schmid value
 twin{tnum}.eta1=Miller(3,-1,2,twin{tnum}.CS,'uvw'); %What you specify here affects sign and schmid value
 twin{tnum}.actType=2; %i.e. 180 around K1, 180 around eta, compound K1 then eta1, compound eta1 then K1
 twin{tnum}.variantsToUse=1; 
 twin{tnum}.tol.misGb=5*degree; %Used to produce the merged grains starting point for segmentation
 twin{tnum}.tol.misGbRlx=10*degree; %Used to 
 twin{tnum}.tol.misMean=10*degree; 
 rot=twin{tnum}.CS.rotation;
 twin{tnum}.variantRotations=rot(:,1);
 twin{tnum}.voteWeights=[2,1,1,1]; %Weights for schmid, boundary, area, and initial ODF
  
 tnum=3;
 twin{tnum}.CS=opt.CS{2};
 twin{tnum}.name='{112}';
 twin{tnum}.k1=Miller(1,1,2,twin{tnum}.CS,'hkl'); %What you specify here affects sign and schmid value
 twin{tnum}.eta1=Miller(3,-7,2,twin{tnum}.CS,'uvw'); %What you specify here affects sign and schmid value
 twin{tnum}.actType=1; %(pick 1-4):180 around K1, 180 around eta1, compound - K1 then eta1, compound - eta1 then K1
 twin{tnum}.variantsToUse=1; 
 twin{tnum}.tol.misGb=5*degree; %Used to produce the merged grains starting point for segmentation
 twin{tnum}.tol.misGbRlx=10*degree; %Used to 
 twin{tnum}.tol.misMean=10*degree; 
 rot=twin{tnum}.CS.rotation;
 twin{tnum}.variantRotations=rot(:,1);
 twin{tnum}.voteWeights=[2,1,1,0]; %Weights for schmid, boundary, area, and initial ODF  based voting
  
 
 %Compute twin properties 
 opt.twin=getTwinProperties(twin);

 %store nTwin and twinUnknown index
 opt.nTwin=length(opt.twin)-1;
 opt.twinUnknown=length(opt.twin);

  %Set Family reconstruction parameter
 opt.grain_recon.FamilyMisTol=12*degree;
 opt.grain_recon.segAngle=6*degree;
 opt.grain_recon.segType=length(opt.twin)+1;

 %Set global merging parameters
 %Notes: A fragment/cluster is only merged if it is fully internal or
 %satisfies 
 opt.mergeByGeometry.mergeTP=false; %MTEX merges all grains at tripple point if two have a twin relationship
 opt.mergeByGeometry.mergeFrag=true; %fragments internal to another fragment are merged
 opt.mergeByGeometry.mergeCluster=true; %clusters or fragments internal to another cluster are merged
 opt.mergeByGeometry.itter=2; %number of itteratation on merging to see if the merged clusters result in more merging cases
 opt.mergeByGeometry.SGT=200; %Small grain threshold for merging consideration using Families in neighboring grains
 opt.mergeByGeometry.ART=3.5; %High aspect ratio merged grain threshold for merging consideration using Families in neighboring grains.
 opt.mergeByGeometry.MARGST=1000; %Largest grain size in pixels for allowing high aspect ratio grainsto be merged
 opt.mergeByGeometry.minSBR=0.7; %Min shared boundary ratio for merging consideration with relaxed twin misorientation tolerance  
 opt.mergeByGeometry.FamilyMisTol=12*degree; %Used to merge grains that satisfy SGT or ART


 %Tree determination 
 %==================================================================================================
 %The root of tree is the parent grain. The first criterion is the centrality metric 
 %outcloseness which is a combined metric of how much of the directional graph can be reached and 
 %the distance to traverse the graph. The directional graph is essentially constructed using the 
 %methodology described in [Marshall 2009]. The second criterion considers whether a grain is in the 
 %initial texture. The weight by volume options uses each fragment fraction in the initial
 %material's ODF to multiple the outcloseness metric. Thus, the parent node must both describe the
 %graph relationships as a whole and be reasonable based on the initial texture. The max iterations
 %will keep the tree from continuing to iterate over the grain and the grain will be returned with
 %an error for the issue. The number of iterations will dependon the number family relationships
 %that have to be determined.
 %==================================================================================================
 opt.rootFamily.weightByFVOL=false;
 opt.familyTree.useMinSpanTree=true;
 opt.familyTree.maxMinSpanTreeItter=40;
 
 
 %Generation weight raises the cost to traverse an edge relationship with a higher generation compared to a lower 
 %generation. This is similar to comparing boundary ratios and schmid in [Marshall 2009]. A value of around 
 %0.5 worked well for a high secondary twinning fraction and small amounts of tertiary twinning. 
 %The formula follows:    edgeweight = (1+genweight*(gen-1))*(edge distance)
 %such that a first generation edge has no modifications and generations 2 has a penalty of
 %1+genweight. This is a weak enforcement that multiple twins of the same mode share the same
 %parent, but is not so stringent that higher order twins of the same mode are completely
 %determined. The same generation has no penalty.
 %==================================================================================================
 opt.familyTree.genWeight=0.5;
 
 %This options enforces that if a Family has twin relationships twins of the same mode but different
 %variants, then that Family is the parent of those twins. This optionswas implmented in ;Marshall
 %2009]. I don't think this is a safe option for when there is third order twinning. 
 opt.familyTree.enforceTwinModeParents=false;  
 
 %Family relationships that don't share boundary can be excluded from consideration. This greatly
 %simplifies the number of circular relationships. Occasionally twin families are removed by
 %enforcing this. To this end there is a list of exception instances where the Family twin
 %relationship should still be included.
 %==================================================================================================
 opt.familyTree.removeNSBE=true; 
 opt.familyTree.addNSBE_NodeRemoved=true;
 opt.familyTree.addNSBE_FavorableParent=true;
 opt.familyTree.noSameType=false;
 
 %The max generation is not enforced explicitly; however, a max edge length is set during the
 %spanning tree when an edge exceeds the generation limit. If the edge is the only explaination for
 %the twin family, then a generation greater than max gen with still occur. Furthermore, a list of grains with
 %generation>maxgen are exported.
 opt.maxGen=3;
 %% Specimen stress state
 opt.sigma=stressTensor([0 0 0; 0 -1 0; 0 0 0]) %Sign of loading does more than just invert twin/parent flag

%% Cleanup local variables 
clear 'tnum' 'twin'
