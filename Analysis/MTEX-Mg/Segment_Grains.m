%% Load the inputs
tic
Inputs


% segment grains
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'angle',5*degree);

% remove two pixel grains
ebsd(grains(grains.grainSize<=2)) = [];
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'angle',5*degree);

% smooth them
grains = grains.smooth(5);

time.inputs=toc;
%% Cluster grains based on fragment misorientation and compute grain fragment properties
% Use a graph to handle and visualize the data. The motivation for using a
% graph is data management and visualization.
tic
doplot=true;

G_Initial = InitializeGraph(ebsd,grains,twin,Mistol,meanMistol,meanMistolRelaxed,doplot);
time.InitializeGraph=toc;
%% Create user edgeList to pre-filter twin relationships
tic
doplot=true;
dolabel=true;
doTwinGrainOnly=true;
toremove=load('toremove.txt')
toadd=load('toadd.txt')
G_clust = ClusterGrainsTwins(G_Initial,grains,Mistol,meanMistolRelaxed,...
    twin,toremove,toadd,doplot,dolabel,doTwinGrainOnly);
time.ClusterGrainsTwins=toc;
%% Compute Schmid info for twin/parents in clustered grains
%This computes Schmid factor for twin/parent identification
tic
G_clust = GetSchmidRelative(G_clust,twin,sigma);
time.GetSchmidRelative=toc;
%% Identify Families in grain clusters 
%Remove nodes that aren't connected by edges, alternatively use condition
%on minimum occurance from conncomp
tic
doplot=true;
dolabel=true;
doTwinGrainOnly=false;
G_clust = AssignFamilyIDs(G_clust,grains,seg_angle_grouped,doplot,dolabel,doTwinGrainOnly);
time.AssignFamilyIDs=toc;
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
tic
w=[1,1,1]
G_clust = FamilyVotes(G_clust,grains,w);
time.FamilyVotes=toc;
%% Prefilter
%% Cleanup the family tree
%The section is the main cleanup of incorrectly characterized twin
%relationships.

% Twin relationship filters
%Apply filters to twin modes
tFilter.twinModes=[1,2,3];
tFilter.singleFragRelationship=true; %only consider removing if twin relationship is with one fragment
tFilter.interalMerge=true; %If grain does not have twin relationship but is internal to a grain set label twin as -1

%Probability of twin variant occuring according schmid variant ranking
tFilter.Schmid.use=false;
tFilter.Schmid.value=6; %not a twin if rank is less than or equal to value
tFilter.Schmid.onlyOutside=true; %only consider removing if twin relationship is not a grain embedded in another grain

%Remove twins that significantly increase the Percentile Average Relative 
%Indented Surface (i.e. if removing a twin improves convexity of the grain cluster
%by some percentage then it should be removed.)
tFilter.PARIS.use=false;
tFilter.PARIS.value=10; %Percent difference improvement in convexity
tFilter.PARIS.onlyOutside=true; %only consider removing if twin relationship is not a grain embedded in another grain

%Remove twin relationships where the ratio of twin boundary of fragment to
%total boundary is less than some minimum ratio
tFilter.RB.use=false; 
tFilter.RB.value=0.02;
tFilter.RB.onlyOutside=true; %only consider removing if twin relationship is not a grain embedded in another grain

%Max grain size twin filter which is taken as a grain size n standard 
%deviations from mean 
tFilter.MaxGS.use=false; 
tFilter.MaxGS.value=4; %n standard deviations above the mean. Not neccisarily an integer
tFilter.MaxGS.onlyOutside=true; %only consider removing if twin relationship is not a grain embedded in another grain

%Min aspect ratio
tFilter.AR.use=false;
tFilter.AR.value=6; %n standard deviations above the mean. Not neccisarily an integer
tFilter.AR.onlyOutside=true; %only consider removing if twin relationship is not a grain embedded in another grain

%Assign unknown twin 
tFilter.ForceInternalunknownTwin.use=false; 
G_Complete_unclean=G_clust;
[G_Complete_unclean]= Filter(G_Complete_unclean,grains,tFilter)


%%
runCleanup=true
maxIter=5;
CleanupIter=0;
%I think this shouldn't take more than 2 itterations. There can be some
%cases that aren't solved until grains are regroupd (i.e you have two
%seperate grains that need still need to be split and a circular relationship exists
%because the same families are in both grains. It's a thing.. go figure!)
time.CleanFamilyTree=0;
G_Complete_unclean=G_clust;
while runCleanup && CleanupIter<maxIter
    tic
    [G_Complete_unclean,runCleanup] = CleanFamilyTree(G_Complete_unclean,grains,tFilter);
    time.CleanFamilyTree=time.CleanFamilyTree+toc;
    if runCleanup
        tic
        G_Complete_unclean=ClusterGrainsTwins(G_Complete_unclean,grains,Mistol,meanMistolRelaxed,...
        twin,[],[],false,false,false);
        time.ClusterGrainsTwins=time.ClusterGrainsTwins+toc;
    end
    
    CleanupIter=CleanupIter+1;
end
if CleanupIter == maxIter
    disp('More CleanupIter than expected, likely need to solve something manually or improve the code')
end
%% Reconstruct grains based on cleanup (some new grains may have been created
 % number of edge from label 
tic
toremove=load('toremove2.txt')
toadd=load('toadd2.txt')
doplot=true;
dolabel=true;
doTwinGrainOnly=true;
G_Complete = ClusterGrainsTwins(G_Complete_unclean,grains,Mistol,meanMistolRelaxed,...
    twin,toremove,toadd,doplot,dolabel,doTwinGrainOnly)
time.ClusterGrainsTwins=time.ClusterGrainsTwins+toc;
%% Redo Families in grain clusters 
%Remove nodes that aren't connected by edges, alternatively use condition
%on minimum occurance from conncomp
doplot=true;
dolabel=false;
doTwinGrainOnly=true;
G_Complete = AssignFamilyIDs(G_Complete,grains,seg_angle_grouped,doplot,dolabel,doTwinGrainOnly);
time.AssignFamilyIDs=time.AssignFamilyIDs+toc;
%% Create the family tree 
%uses the recursive relationship as described in Pradalier et al. 2018
tic
[G_Complete] = CreateFamilyTree(G_Complete,grains);
time.CreateFamilyTree=toc;
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
[grains_fin] = transferGtoGrains(G_Complete2,grains)
save('grains_fin.mat', 'grains_fin');
time.transferGtoGrains=toc;
%% Final plots
figure; plot(grains_fin,grains_fin.prop.Type)
print('twin_type','-dtiffn','-r300');
figure; plot(grains_fin,grains_fin.prop.Generation)
print('twin_gen','-dtiffn','-r300');
figure; histogram(grains_fin.prop.schmidActiveRank(grains_fin.prop.schmidActiveRank>0),6,'BinLimits',[0.5,6.5],'BinMethod','integers');
xlabel('Variant schmid rank');ylabel('Count')
print('twin_schmid_rank','-dtiffn','-r300');
%% Timing 
time
timetotal=struct2cell(time);
timetotal=sum([timetotal{:}])/60
