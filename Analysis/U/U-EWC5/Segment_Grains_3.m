%% Load the inputs
Inputs;
load('EBSDx3_Clean.mat');
load('Grainsx3_Clean.mat');
%% Build master graph and compute preliminary clusters based on boundary misorientation
% Use a graph to handle and visualize the data. 
opt.plot.do=true;
opt.plot.ClusterOnly=true;
opt.plot.labelEdges=false;
opt.plot.legendOn=false;
[G,mGbGrains,mGbGrainsRlx] = InitializeGraph(ebsd,grains,opt);
figure;plot(ebsd,ebsd.orientations)

%% Remove false twin edges and add whatever is missing

opt.plot.do=false;
opt.plot.ClusterOnly=true;
opt.plot.labelEdges=true;
opt.plot.legendOn=false;
[G_clust,G,mGrains] = Cluster(G,grains,opt);

figure;plot(grains,grains.meanOrientation);hold on;
plot(mGrains.boundary,'lineWidth',2,'lineColor','k');hold off

figure;plot(grains,G_clust.Nodes.FamilyID);hold on;
plot(mGrains.boundary,'lineWidth',2,'lineColor','w');hold off
%% Produces a list of groups with no aparent parent also build the family relation matrix.
%This function should return clean before going on to CreateFamilyTree,
%though one can proceed and CreateFamilyTree may be able to handle input by
%setting the grain to unsolved. 
[G_clust,time]= CleanFamilyTree(G_clust,grains,mergedGrains,twin,seg_angle_grouped,true,3,time);
%% Create the family tree 
%uses the recursive relationship as described in Pradalier et al. 2018
tic
[G_clust,err_group] = CreateFamilyTree(G_clust,grains,mergedGrains,twin)
time.CreateFamilyTree=toc;
figure; plot(grains,G_clust.Nodes.Type);mtexColorbar;mtexTitle('Twin Type')
hold on;plot(mergedGrains.boundary,'linecolor','w','linewidth',2.5,'linestyle','-','displayName','merged grains');hold off
text(mergedGrains,int2str(mergedGrains.id));hold off;

figure; plot(grains,G_clust.Nodes.Generation); mtexColorbar;;mtexTitle('Generation')
hold on;plot(mergedGrains.boundary,'linecolor','w','linewidth',2.5,'linestyle','-','displayName','merged grains');hold off

% text(mergedGrains,int2str(mergedGrains.id));hold off;
%% Look at grain clusters with large number of generations
% groups=unique([])
% value=G_Complete.Nodes.Type;
% plotNeighbors=true;
% enforceClusterOnlyMod=false;
% [G_clust,~,~] = ClusterEditor(groups,G_clust,grains,mergedGrains,value,0,plotNeighbors,enforceClusterOnlyMod); 
%% Save the graph and input variables
G_clust.Nodes.FgB=[]; %Can this be removed?
save('Segmentation.mat', 'G','G_clust','ebsd','grains','mergedGrains','CS',...
    'meanMistol','meanMistolRelaxed','Mistol','MistolRlx','min_CI',...
    'min_grainSz','minNEdgeMistol','seg_angle','seg_angle_grouped',...
    'sigma','twin','voteWeights');

