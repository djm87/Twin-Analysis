function [G_clust,mergedGrains,time] = ClusterGrainsTwins(G,grains,type,typeRlx,...
    meanMistol,meanMistolRelaxed,minNEdgeMistol,twin,doPlot,...
    plotClusterOnly,time)
    %To add or remove edge specify the id in eAddList.txt or eRemoveList.txt
    %To add a node and try all edges with relaxed twin tolerance specify node
    %id in nAddList.txt
    %To remove all edges associated with a node, specify node id in
    %nRemoveList.txt
    %Edges are automatically added to grain internal boundaries and twin
    %relationships are tested.
    %All other edge clusters are determined by MergeByBoundary and filtered
    %using the relaxed mean orientation relationships.
    
    tic
    
    eRemoveList=load('eRemoveList.txt'); %Remove the edge 
    eAddList=load('eAddList.txt'); %Adds edge and tries to give it a twin relationship, otherwise type twin unknown
    nRemoveList=load('nRemoveList.txt'); %Removes all edges associated with node
    nAddList=load('nAddList.txt'); %Tries to add all edges that have twin relationship
        
    G.Edges.type=typeRlx;
    
    %Any edges that aren't a twin should be set to be removed
    %Filter edge by mistol based with the type of meanMistolRelaxed
    toRemove=~G.Edges.combineBoundary;
    toRemove(typeRlx==0)=true;
    
    %Initialize combine boundary from boundary merge
    G.Edges.combineCleaned = UseMeanForSmallBoundary(G.Edges.pairs,G.Edges.combineBoundary,grains,type,minNEdgeMistol);
    G.Edges.combineCleaned(toRemove)=false; %
    
    %Try adding nodes
    for i=1:length(nAddList)
        [~,epairs]=grains(nAddList(i)).neighbors;
        for j=1:size(epairs,1)
            eId=find(all(epairs(j,:)==G.Edges.pairs,2) | all(fliplr(epairs(j,:))==G.Edges.pairs,2));
            isTwin=G.Edges.type(eId)~=0 & G.Edges.type(eId)~=length(twin);
            G.Edges.combineCleaned(eId(isTwin))=true;
        end
    end
    
    
    %Try removing nodes
    for i=1:length(nRemoveList)
        [~,epairs]=grains(nRemoveList(i)).neighbors;
        for j=1:size(epairs,1)
            eId=find(all(epairs(j,:)==G.Edges.pairs,2) | all(fliplr(epairs(j,:))==G.Edges.pairs,2));
            isTwin=G.Edges.type(eId)~=0 & G.Edges.type(eId)~=length(twin);
            G.Edges.combineCleaned(eId(isTwin))=false;
        end
    end
    
    %Add edges
    typeAdd=typeRlx(eAddList);
    typeAdd(typeAdd==0)=length(twin); %Unknown Type
    G.Edges.type(eAddList)=typeAdd;
    G.Edges.combineCleaned(eAddList)=true;
    
    %Update unknown types in nodes
    typeKnown=unique(G.Edges.pairs([G.Edges.type~=length(twin),G.Edges.type~=length(twin)]));
    G.Nodes.typeUnknown(:)=true;
    G.Nodes.typeUnknown(typeKnown)=false;
    
    %Remove edges
    G.Edges.combineCleaned(eRemoveList)=false; %
    
    %Merge grains based on edge combine list
    [mergedGrains,parentId,mergedTwinBoundary,G.Edges.combineCleaned] = MergeByEdge(G.Edges.pairs,G.Edges.combineCleaned,grains);    

    %Remove edges to make reduced graph over the clusters
    G_clust=rmedge(G,G.Edges.pairs(~G.Edges.combineCleaned,1),...
        G.Edges.pairs(~G.Edges.combineCleaned,2));
    
    %Assign group Ids
    G_clust = GroupsFromGraph(G_clust,mergedGrains,parentId);
    
    %Plot the results
    if doPlot==true
        %Remove grains that haven't twinned (better plot performance)
        if plotClusterOnly
            tokeep=zeros(length(G_clust.Nodes.Id),1,'logical');
            tokeep(unique(G_clust.Edges.pairs))=true;
            G_clust_clean=rmnode(G_clust,find(~tokeep));
        else
            G_clust_clean=G_clust;
        end
        
        %Plot edge labeled graph
        labelNodes=false;labelEdges=true;plotG=true;legendOn=false;
        fhandle = plotGraph(grains,mergedGrains,G_clust_clean,...
            grains.meanOrientation,G_clust_clean.Nodes.Id,...
            labelNodes,labelEdges,legendOn,plotG,[]);
        
        %Plot node labeled graph
        labelNodes=true;labelEdges=false;plotG=false;legendOn=false;
        fhandle = plotGraph(grains,mergedGrains,G_clust_clean,...
            grains.meanOrientation,G_clust_clean.Nodes.Id,...
            labelNodes,labelEdges,legendOn,plotG,[]);
        
    end
    if ~isfield(time,'ClusterGrainsTwins')
        time.ClusterGrainsTwins=0;
    end
    time.ClusterGrainsTwins=time.ClusterGrainsTwins+toc;
    
end

