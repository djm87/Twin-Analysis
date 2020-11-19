%% Load the inputs
Inputs;
load('EBSDx3_Clean.mat');
load('Grainsx3_Clean.mat');

%% Cluster grains based on fragment misorientation and compute grain fragment properties
% Use a graph to handle and visualize the data. 
doPlot=false;
doLabel=false;
[G_Initial,time] = InitializeGraph(ebsd,grains,twin,Mistol,meanMistol,...
    meanMistolRelaxed,doPlot,doLabel,time);

%% Do standard MTEX clustering and mean orientation
%The use of mean orientations require a relaxed tolerance for twin 
%relationship identification and often miss relationships if the
%misorientation in a grain is significant. However, grain to grain
%relationships are what we are interested in, not boundaries. To address
%this, the boundary merging method is used for initial clustering, while
%the type of the relationship is determined by significantly relaxing the 
%mean orientation relationships.
doPlot=true;
[G_Initial.Edges.combineBoundary,mergedGrains,twinBoundary,time] = ...
    MergeByBoundary(G_Initial,grains,Mistol,twin,doPlot,time);
[G_Initial.Edges.combineBoundaryRlx,mergedGrainsRlx,twinBoundaryRlx,time] = ...
    MergeByBoundary(G_Initial,grains,MistolRlx,twin,doPlot,time);


%Get the type of misorientation
mori=inv(grains(G_Initial.Edges.pairs(:,1)).meanOrientation).*...
    grains(G_Initial.Edges.pairs(:,2)).meanOrientation; 
[~,tType] = TestTwinRelationship(mori,meanMistol,twin,G_Initial.Edges.type);
[~,tTypeRlx] = TestTwinRelationship(mori,meanMistolRelaxed,twin,G_Initial.Edges.type);

%% Create user edgeList to pre-filter twin relationships
%To add or remove edge specify the id in eAddList.txt or eRemoveList.txt
%To add a node and try all edges with relaxed twin tolerance specify node
%id in nAddList.txt
%To remove all edges associated with a node, specify node id in
%nRemoveList.txt
%Edges are automatically added to grain internal boundaries and twin
%relationships are tested.
%All other edge clusters are determined by MergeByBoundary and filtered
%using the relaxed mean orientation relationships.

doPlot=false;
plotClusterOnly=false;
[G_clust,mergedGrains,time] = ClusterGrainsTwins(G_Initial,grains,tType,tTypeRlx,...
    meanMistol,meanMistolRelaxed,minNEdgeMistol,twin,doPlot,...
    plotClusterOnly,time);


%% Identify Families in grain clusters 
%Give each cluster a group Id and for each cluster group similar
%orientations into families, giving them an unique Id
doPlot=true;
[G_clust,time]= AssignFamilyIDs(G_clust,grains,mergedGrains,seg_angle_grouped,twin,doPlot,time);

%% Look at grain clusters with a large number of family ids 
groups=unique(G_clust.Nodes.Group(G_clust.Nodes.FamilyID>6));
value=grains.meanOrientation;
plotNeighbors=false;
enforceClusterOnlyMod=false;
% [G_clust,~,~] = ClusterEditor(groups,G_clust,grains,mergedGrains,value,0,plotNeighbors,enforceClusterOnlyMod); 

%% Look at grain clusters with a large number of family ids 
groups=mergedGrains(mergedGrains.aspectRatio>5).id;
value=grains.meanOrientation;
plotNeighbors=true;
enforceClusterOnlyMod=false;
% [G_clust,~,~] = ClusterEditor(groups,G_clust,grains,mergedGrains,value,0,plotNeighbors,enforceClusterOnlyMod); 

%% Look at custom list of group Id 
% figure;plot(grains,G_clust.Nodes.FamilyID);mtexColorMap(hsv);hold on;
% plot(mergedGrains.boundary,'lineWidth',2,'lineColor','k')
% text(mergedGrains,int2str(mergedGrains.id));hold off;
% figure;plot(grains,grains.meanOrientation);hold on;
% plot(mergedGrains.boundary,'lineWidth',2,'lineColor','w')
% text(mergedGrains,int2str(mergedGrains.id));hold off;
%% Breakup
groups=unique([])
value=grains.meanOrientation;
plotNeighbors=true;
enforceClusterOnlyMod=false;
% [G_clust,~,~] = ClusterEditor(groups,G_clust,grains,mergedGrains,value,0,plotNeighbors,enforceClusterOnlyMod); 

%% Find clusters that look like they aren't right

%Get the decrease Percentile Average Relative Indented Surface after removing
% a twin 
% [G_clust,time]=eRemovalPARISChange(G_clust,grains,true,[1,2,3],time)


% groups=unique(G_clust.Nodes.Group(G_clust.Nodes.ePARISChange >.1));
% value=zeros(length(G_clust.Nodes.Id),1);
% value(G_clust.Nodes.ePARISChange >.1)=1;
% plotNeighbors=false;
% enforceClusterOnlyMod=false;
% [G_clust,~,~] = ClusterEditor(groups,G_clust,grains,mergedGrains,value,0,plotNeighbors,enforceClusterOnlyMod); 

%% Compute Schmid info for twin/parents in clustered grains
%This computes Schmid factor for twin/parent identification
doPlot=true;
[G_clust,time]= GetSchmidRelative(G_clust,grains,mergedGrains,twin,sigma,doPlot,time);


%% Implement voting scheme for the respective families
%For each family compute the vote between connected families 
%Edges currently connect clustered grains. We can compute the vote for each
%edge and then average for the family. This ensures that only neighbors are
%ever voted on. So that votes are really between familys average values of
%schmid are computed and stored for each fragment in a family. The area and
%relative boundary are taken in the addative sense rather than average. 

tic
G_clust = FamilyVotes(G_clust,grains,voteWeights);
time.FamilyVotes=toc;

%%
[G_Complete,time]= CleanFamilyTree(G_clust,grains,mergedGrains,twin,seg_angle_grouped,true,3,time);

%% Create the family tree 
%uses the recursive relationship as described in Pradalier et al. 2018
tic
[G_Complete] = CreateFamilyTree(G_Complete,grains,twin);
time.CreateFamilyTree=toc;
figure; plot(grains,G_Complete.Nodes.Type);mtexColorbar;mtexTitle('Twin Type')
hold on;plot(mergedGrains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-','displayName','merged grains');hold off
text(mergedGrains,int2str(mergedGrains.id));hold off;

figure; plot(grains,G_Complete.Nodes.Generation); mtexColorbar;;mtexTitle('Generation')
text(mergedGrains,int2str(mergedGrains.id));hold off;
%% Look at grain clusters with large number of generations
groups=unique(G_Complete.Nodes.Group(G_Complete.Nodes.Generation<0));
value=G_Complete.Nodes.Generation;
plotNeighbors=false;
enforceClusterOnlyMod=false;
% [G_clust,~,~] = ClusterEditor(groups,G_Complete,grains,mergedGrains,value,0,plotNeighbors,enforceClusterOnlyMod); 

%% Extract parent stats for schmid and compute the angle
tic
G_Complete2 = GetSchmidVariants(G_Complete,twin,sigma);
time.GetSchmidVariants=toc;
%% Reconstruct twins to fix count issues
%***Added list to merge for counting
tic
[~,G_Complete2] = fragmentReconstruction(G_Complete2,ebsd,grains);
time.fragmentReconstruction=toc;
figure; plot(grains,G_Complete2.Nodes.MergeTwin)

%% Get count statistics 
tic
[stats_fin.twinCount,G_Complete2] = CountTwins(G_Complete2)
time.CountTwins=toc;
figure;histogram(stats_fin.twinCount(stats_fin.twinCount>0))
xlabel('Number of twins in grain');ylabel('Count')
print('twin_count','-dtiffn','-r300');
% figure; 
% plot(grains,G_Complete2.Nodes.twinCount)
% text(grains,int2str(G_Complete2.Nodes.twinCount))

%% Get twin thickness 
tic
[G_Complete2] = TwinThickness(G_Complete2,grains,twin)
time.TwinThickness=toc;
% selection=G_Complete2.Nodes.twinThickness<1;
% figure; [~,mP] = plot(grains(selection),G_Complete2.Nodes.twinThickness(selection))
% mP.micronBar.length=2;

%% Compute the volume fractions 
tic
[stats_fin.twinVF,stats_fin.totalArea,stats_fin.totalTwinVF] = ...
    getTwinFractions(G_Complete2,grains,twin)
save('stats_fin.mat', 'stats_fin');
time.getTwinFractions=toc;
%% Transfer G to grains
tic
[grains_fin,mergedGrains_fin] = transferGtoGrains(G_Complete2,grains,mergedGrains,twin)
save('grains_fin.mat', 'grains_fin');
time.transferGtoGrains=toc;
%% Final plots
figure; plot(grains_fin,grains_fin.prop.Type);hold on;
plot(mergedGrains.boundary,'lineWidth',2,'lineColor','w')
text(mergedGrains,int2str(mergedGrains.id));hold off;
print('twin_type','-dtiffn','-r300');
figure; plot(grains_fin,grains_fin.prop.Generation)
print('twin_gen','-dtiffn','-r300');

AW=grains_fin.area./stats_fin.totalArea;

uniqueGen=unique(grains_fin.prop.Generation(grains_fin.prop.Generation>=0));
genArea=zeros(length(uniqueGen),1);
for i=1:length(uniqueGen)
    rankAW=AW(grains_fin.prop.Generation==uniqueGen(i));
    genArea(i)=sum(rankAW);
end
figure; bar(uniqueGen,genArea)

AW=grains_fin.area./stats_fin.totalArea;
rankBar=zeros(6,1);
for i=1:6
    rank=grains_fin.prop.schmidActiveRank(grains_fin.prop.schmidActiveRank==i);
    rankAW=AW(grains_fin.prop.schmidActiveRank==i);
    rankBar(i)=sum(rankAW);
end
figure;bar(rankBar)

groups=grains_fin.prop.Group;
minRankGroup=zeros(length(groups),1);
nRankGroup=zeros(length(groups),1);
for i=1:length(groups)
    group=groups(i);
    nId=find(group==grains_fin.prop.Group);
    schmidActiveRank=grains_fin.prop.schmidActiveRank(nId);
    schmidActiveRankNZ=schmidActiveRank(schmidActiveRank>0);
    uniqueRank=unique(schmidActiveRankNZ);
    if uniqueRank~=0
        nRankGroup(i)=length(uniqueRank);
        minRankGroup(i)=min(schmidActiveRankNZ);
    end
end
figure; histogram(nRankGroup,6,'BinLimits',[0.5,6.5],'BinMethod','integers');
figure; histogram(minRankGroup,6,'BinLimits',[0.5,6.5],'BinMethod','integers');
figure; scatter(nRankGroup,grains_fin.area)
figure; scatter(minRankGroup,grains_fin.area)
figure; scatter(minRankGroup,mergedGrains(grains_fin.prop.Group).area)

indx=grains_fin.prop.Generation==1;
xval=grains_fin.prop.totalForTwinnedArea(indx)./grains_fin.prop.totalForAllArea(indx);
figure; scatter(mergedGrains_fin.prop.twinArea./mergedGrains_fin.area,mergedGrains_fin.diameter)
figure; histogram(grains_fin.prop.schmidActiveRank(grains_fin.prop.schmidActiveRank>0),6,'BinLimits',[0.5,6.5],'BinMethod','integers');
xlabel('Variant schmid rank');ylabel('Count')
print('twin_schmid_rank','-dtiffn','-r300');
%% Timing 
time
timetotal=struct2cell(time);
timetotal=sum([timetotal{:}])/60
