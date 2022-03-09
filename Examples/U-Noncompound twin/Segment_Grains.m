%% Load the inputs/data

Inputs;

%% Build inital fragment graph 
[G,mGbGrains] = InitialGraph(ebsd,grains,opt);

%Plot full graph
labelNodes=true;labelEdges=false;plotG=true;legendOn=false;
fhandle = plotGraph(grains,[],G,...
    grains.meanOrientation,G.Nodes.Id,...
    labelNodes,labelEdges,legendOn,plotG,[]);
mtexTitle('Full graph')      

%Plot boundary based clusters
figure;plot(grains,grains.meanOrientation,'noBoundary');hold on;
plot(mGbGrains.boundary,'lineWidth',2,'lineColor','w');
% text(mGbGrains,int2str(mGbGrains.id))
mtexTitle('Boundary based merged grains');hold off  


%% Build cluster graph by removing, adding, etc... edges

[G_clust,G,mGrains] = ClusterGraph(G,grains,opt);

%Plot mGrains
figure;plot(grains,grains.meanOrientation,'noBoundary');hold on;
plot(mGrains.boundary,'lineWidth',2,'lineColor','w');
text(mGrains,int2str(mGrains.id))
mtexTitle('Cluster graph merged grains');hold off 

%% Edit merged clusters

%plot quantities such as FamilyID that can identify problem clusters
figure;plot(grains,G.Nodes.FamilyID,'noBoundary');hold on;
plot(mGrains.boundary,'lineWidth',2,'lineColor','w');
text(mGrains,int2str(mGrains.id));hold off;
mtexTitle('Cluster graph Families')  


% groups=unique(G_clust.Nodes.Group(G_clust.Nodes.FamilyID>5));
% groups=[115 122 80]
groups=[]
value=G.Nodes.FamilyID;
GraphEditor(groups,1,[],G_clust,G,grains,mGrains,value,0,1,0,1,2);


%% Build family graph for each cluster

computemGrainId=[];
% computemGrainId=mGrains.id;
[G_Family,G_clust,G] = FamilyGraph(G_clust,G,grains,computemGrainId,opt);


%% Clean the family tree

cleanGroups=unique(G_Family.Nodes.Group);
[G_Family,G_clust,G,exflagGroup]= CleanFamilyTree(cleanGroups,...
    G_Family,G_clust,G,grains,mGrains,opt);   


%% Visualize the results

%Plot grain fragment twin type
figure; plot(grains,G.Nodes.type,'noBoundary');
hold on;plot(mGrains.boundary,'linecolor','w','linewidth',2.5,'linestyle','-','displayName',...
    'merged grains');hold off
mtexColorbar;mtexTitle('Twin Type');
saveFigure('Twin type')

%Plot the number of generations in each cluster
figure; plot(grains,G.Nodes.Generation,'noBoundary');
hold on;plot(mGrains.boundary,'linecolor','w','linewidth',2.5,'linestyle','-','displayName',...
    'merged grains');hold off;
mtexColorbar;mtexTitle('Generation');
saveFigure('Twin Generation')

figure; plot(grains,G.Nodes.EffSF,'noBoundary');
hold on;plot(mGrains.boundary,'linecolor','w','linewidth',2.5,'linestyle','-','displayName',...
    'merged grains');hold off;
mtexColorbar;mtexTitle('EFF SF');

figure; plot(grains,G.Nodes.nSFAV,'noBoundary');
hold on;plot(mGrains.boundary,'linecolor','w','linewidth',2.5,'linestyle','-','displayName',...
    'merged grains');hold off;
mtexColorbar;mtexTitle('SF active variant');

figure; plot(grains,G.Nodes.nSFAVR,'noBoundary');
hold on;plot(mGrains.boundary,'linecolor','w','linewidth',2.5,'linestyle','-','displayName',...
    'merged grains');hold off;
mtexColorbar;mtexTitle('Active variant SF rank (lower is most stresses)');

%% Edit/Visualize the Family tree

%Use the exit flag to construct groups to manually address
%exflagGroup=-1 single frag
%exflagGroup=1 no twin relationships but multiple families
%exflagGroup=2 too many parents
%exflagGroup=3 no parent but twin relationships exist
%exflagGroup=4 max number of generations hit while assigning generation
%exflagGroup=5 not all fragments were related... something is wrong
%exflagGroup=6 entered fix circular relationship fnc too many times
groups=cleanGroups(exflagGroup>0)
% groups=unique(G_clust.Nodes.Group(find(G_clust.Nodes.computeFamily)))
groups=[]
value=G.Nodes.FamilyID;
[G_Family] = GraphEditor(groups,1,G_Family,G_clust,grains,mGrains,value,0,1,0,1,1);


%% Compute the twin fraction 

[mGrains,twinVF] = getTwinFractions(G,grains,mGrains,opt);


%% Twin thickness 

[G] = TwinThickness(G,grains,opt);
figure;histogram(G.Nodes.twinThickness(G.Nodes.twinThickness>0));
title('# Twin Fragment thickness');



%% Twin count

[mGrains] = CountTwins(G,grains,mGrains,opt);

figure;histogram(mGrains.prop.twinCount,'BinMethod','integers');title('# Twin Fragments per cluster');
figure;histogram(mGrains.prop.twinFamilyCount,'BinMethod','integers');title('# Unique Twin Families per cluster');xticks([0,1,2]); 


%% Twin variant ranking

figure;histogram(G.Nodes.nSFAVR(G.Nodes.nSFAVR>0),'BinLimits',[0.5,6.5],'BinMethod','integers');title('# Twin Schmid Rank');xticks([1,2,3,4,5,6]); 


%% Transfer results to grains 

[grains] = transferGtoGrains(G,grains);


%% Save important variable
%save('Segmented_data.mat','G','grains','mGrains','opt','twinVF');
