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
[G,mGbGrains,mGbGrainsRlx] = InitializeGraph(ebsd,grains,opt);
% profile viewer
% profile off
figure;plot(grains,grains.meanOrientation);hold on;
plot(mGbGrains.boundary,'lineWidth',2,'lineColor','w');hold off
%% Remove false twin edges and add whatever is missing
load('cluster.mat')

% profile on
opt.plot.do=false;
opt.plot.ClusterOnly=true;
opt.plot.labelEdges=false;
opt.plot.legendOn=false;
[G_clust,G,mGrains] = Cluster(G,grains,opt);
% profile viewer
% profile off
% profile on
%%

figure;plot(grains,grains.meanOrientation,'noBoundary');hold on;
% figure;plot(ebsd,ebsd.orientations);hold on
plot(mGrains.boundary,'lineWidth',2,'lineColor','w');hold off
text(mGrains,int2str(mGrains.id));hold off;
% profile viewer
% profile off

% save('cluster.mat','G_clust','G','mGrains','mGbGrains','mGbGrainsRlx','grains')
%% Prefilter by addressing problem clusters
groups=unique([1284 1272 1261 1258 1281 1295 1278 1298 1316 1277 1248 1304 1299 1239 1254 82 1323 1339 1372 261 1375 1429 1386 1379 1404 1528 1536 1500 1501 1451 1445 1544 1553 1723 1642 1635 1601 1576 1543 1610 1586 1569 806 1596 1695 1733 1737 1697 1719 1796 1671]);
groups=groups(groups>1418);
groups=[1273];

% figure;plot(grains,G_clust.Nodes.FamilyID)
groups=unique(G_clust.Nodes.Group(G_clust.Nodes.FamilyID>10));


groups=unique([1262 1242 1201 1224 1266 1252 1287 1227 1221 1210 1175 1218 1233 1245 1279 1238 1217 1187 1255 1216 1226 1256 1161 1211 1177 1259 1273 1289 1315 1298]);
groups=groups(groups>=1262);
length(groups)
value=grains.meanOrientation;
[G_clust,~,~] = ClusterEditor(groups,G_clust,grains,mGrains,value,0,1,1,0,1,0); 

%% Produces a list of groups with no aparent parent also build the family relation matrix.
%This function should return clean before going on to CreateFamilyTree,
%though one can proceed and CreateFamilyTree may be able to handle input by
%setting the grain to unsolved. 
CleanFamilyTree(G_clust,grains,mGrains,opt)   

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
groups=unique([1070])
value=grains.meanOrientation;
plotNeighbors=true;
enforceClusterOnlyMod=false;
[G_clust,~,~] = ClusterEditor(groups,G_clust,grains,mGrains,value,0,plotNeighbors,enforceClusterOnlyMod); 
%% Save the graph and input variables
G_clust.Nodes.FgB=[]; %Can this be removed?
save('Segmentation.mat', 'G','G_clust','ebsd','grains','mergedGrains','CS',...
    'meanMistol','meanMistolRelaxed','Mistol','MistolRlx','min_CI',...
    'min_grainSz','minNEdgeMistol','seg_angle','seg_angle_grouped',...
    'sigma','twin','voteWeights');

