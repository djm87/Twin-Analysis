%% Load the inputs
Inputs;
load('EBSD_Clean.mat');
load('Grains_Clean.mat');
%% Build master graph and compute preliminary clusters based on boundary misorientation
% Use a graph to handle and visualize the data. 
% profile on
opt.plot.do=true;
opt.plot.ClusterOnly=false;
opt.plot.labelEdges=false;
opt.plot.legendOn=false;
[G,mGbGrains,mGbGrainsRlx] = InitialGraph(ebsd,grains,opt);
% profile viewer
% profile off
figure;plot(grains,grains.meanOrientation);hold on;
plot(mGbGrains.boundary,'lineWidth',2,'lineColor','w');hold off
%% Build final cluster graphs by removing, adding, etc... edges
% profile on
opt.plot.do=false;
opt.plot.ClusterOnly=false;
opt.plot.labelEdges=false;
opt.plot.legendOn=false;
[G_clust,G,mGrains] = ClusterGraph(G,grains,opt);
% profile viewer
% profile off
% profile on

%% Edit merged clusters
% figure;plot(grains,G_clust.Nodes.FamilyID)
% groups=unique(G_clust.Nodes.Group(G_clust.Nodes.FamilyID>10));

groups=[];
% % groups=groups(groups>=451);
% length(groups)
% value=grains.meanOrientation;
[G_clust,~,~] = ClusterEditor(groups,G_clust,grains,mGrains,value,0,1,1,0,1,0); 
 

%% Build family graph for each cluster
computemGrainId=[];
[G_Family,G_clust,G] = FamilyGraph(G_clust,G,grains,computemGrainId,opt)









%% Build the family tree
%Move to cluster or new graph function
opt.debugFamilyTree=false;
groups=mGrains.id;
exflagGroup=zeros(max(mGrains.id),1);
[G_Family,G_clust,exflagGroup]=CleanFamilyTree(groups,...
    G_Family,G_clust,G,grains,mGrains,exflagGroup,opt);

%% Debug the family tree
opt.debugFamilyTree=false;
groups=find(exflagGroup==2);
exflagGroup=zeros(max(mGrains.id),1);
[G_Family,G_clust,exflagGroup]=CleanFamilyTree(groups,...
    G_Family,G_clust,G,grains,mGrains,exflagGroup,opt);


%%


figure; plot(grains,G_clust.Nodes.Type,'noBoundary');mtexColorbar;mtexTitle('Twin Type')
hold on;plot(mGrains.boundary,'linecolor','w','linewidth',2.5,'linestyle','-','displayName','merged grains');hold off
hold off;

figure; plot(grains,G_clust.Nodes.Generation,'noBoundary');mtexColorbar;mtexTitle('Generation')
hold on;plot(mGrains.boundary,'linecolor','w','linewidth',2.5,'linestyle','-','displayName','merged grains');hold off
hold off;





%%
%%
%%
list=mGrains.grainSize>200

sum(mGrains.diameter.*mGrains.area/sum(mGrains.area))

% groups=unique(G_clust.Nodes.Group(G_clust.Nodes.>10));

figure;plot(grains,grains.meanOrientation,'noBoundary');hold on;
% figure;plot(ebsd,ebsd.orientations);hold on
plot(mGrains.boundary,'lineWidth',2,'lineColor','w');hold off
text(mGrains,int2str(mGrains.id));hold off;

figure;plot(grains,G_clust.Nodes.FamilyID,'noBoundary');hold on;
% figure;plot(ebsd,ebsd.orientations);hold on
plot(mGrains(list).boundary,'lineWidth',2,'lineColor','w');hold off
text(mGrains(list),int2str(mGrains(list).id));hold off;
% profile viewer
% profile off

% save('cluster.mat','G_clust','G','mGrains','mGbGrains','mGbGrainsRlx','grains')


%% Remove the grains along the boarder of the map 
% this gives the boundary ids that are cut
boundaryIds = mGrains.boundary.hasPhaseId(0);

% the corresponding grain ids are
grainId = unique(mGrains.boundary(boundaryIds).grainId(:,2))

% remove those grains from the list
mGrains2=mGrains;
mGrains2(grainId) = []


figure;plot(grains,grains.meanOrientation,'noBoundary');hold on;
plot(mGrains2.boundary,'lineWidth',2,'lineColor','w');hold off
text(mGrains2,int2str(mGrains2.id));hold off;
%% Produces a list of groups with no aparent parent also build the family relation matrix.
%This function should return clean before going on to CreateFamilyTree,
%though one can proceed and CreateFamilyTree may be able to handle input by
%setting the grain to unsolved. 



find(exflagGroup==6)
length(problemClusters)

%% Create the family tree 
%uses the recursive relationship as described in Pradalier et al. 2018
[G_clust2,err_group] = CreateFamilyTree(G_clust,grains,mGrains2,opt)

figure; plot(grains,G_clust.Nodes.Type,'noBoundary');mtexColorbar;mtexTitle('Twin Type')
hold on;plot(mGrains.boundary,'linecolor','w','linewidth',2.5,'linestyle','-','displayName','merged grains');hold off
hold off;
text(mergedGrains,int2str(mergedGrains.id));hold off;

figure; plot(grains,G_clust2.Nodes.Generation); mtexColorbar;;mtexTitle('Generation')
hold on;plot(mergedGrains.boundary,'linecolor','w','linewidth',2.5,'linestyle','-','displayName','merged grains');hold off

% text(mergedGrains,int2str(mergedGrains.id));hold off;
%% Save the graph and input variables
G_clust.Nodes.FgB=[]; %Can this be removed?
save('Segmentation.mat', 'G','G_clust','ebsd','grains','mergedGrains','CS',...
    'meanMistol','meanMistolRelaxed','Mistol','MistolRlx','min_CI',...
    'min_grainSz','minNEdgeMistol','seg_angle','seg_angle_grouped',...
    'sigma','twin','voteWeights');

