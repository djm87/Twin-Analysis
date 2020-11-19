%% Load the inputs
load('Segmentation.mat');

%% Extract parent stats for schmid and compute the angle
tic
G_clust = GetSchmidVariants(G_clust,twin,sigma);
time.GetSchmidVariants=toc;
%% Reconstruct twins to fix count issues
%***Added list to merge for counting
tic
[~,G_clust] = fragmentReconstruction(G_clust,ebsd,grains);
time.fragmentReconstruction=toc;
figure; plot(grains,G_clust.Nodes.MergeTwin)

%% Get count statistics 
tic
[stats.twinCount,G_clust] = CountTwins(G_clust)
time.CountTwins=toc;
figure;histogram(stats.twinCount(stats.twinCount>0))
xlabel('Number of twins in grain');ylabel('Count')
print('twin_count','-dtiffn','-r300');
% figure; 
% plot(grains,G_Complete.Nodes.twinCount)
% text(grains,int2str(G_Complete.Nodes.twinCount))

%% Get twin thickness 
tic
[G_clust] = TwinThickness(G_clust,grains,twin)
time.TwinThickness=toc;
% selection=G_Complete.Nodes.twinThickness<1;
% figure; [~,mP] = plot(grains(selection),G_Complete.Nodes.twinThickness(selection))
% mP.micronBar.length=2;

%% Compute the volume fractions 
tic
[stats.twinVF,stats.totalArea,stats.totalTwinVF] = ...
    getTwinFractions(G_clust,grains,twin)
save('stats.mat', 'stats');
time.getTwinFractions=toc;
%% Transfer G to grains
tic
[grains,mergedGrains] = transferGtoGrains(G_clust,grains,mergedGrains,twin)
time.transferGtoGrains=toc;
save('Result.mat','ebsd','grains','mergedGrains','CS','stats');

