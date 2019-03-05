function G_clust = ClusterGrainsTwins(G,grains,Mistol,meanMistolRelaxed,twin,...
    removeList,doplot,dolabel)
    
%Clusters grains by removing edges of grains that should not be combined
    %Returns G.Edges.combineBoundary and G.Edges.groupMerge (currently not
    %returned)
    [G,mergedGrains]= MergeByBoundary(G,grains,Mistol,twin) 
    
    %Compute misorientation between all pairs
    mori=inv(grains(G.Edges.pairs(:,1)).meanOrientation).*...
        grains(G.Edges.pairs(:,2)).meanOrientation; 

    %Get the type of misorientation
    [~,type] = TestTwinRelationship(mori,meanMistolRelaxed,twin);
    G.Edges.type=type;

    %Any grains that don't match up with meanMistolRelaxed should be
    %removed
    toRemove=G.Edges.type==0;      
    
    G.Edges.combineCleaned=G.Edges.combineBoundary; %Reinitialize combine
    G.Edges.combineCleaned(removeList)=false; %
    G.Edges.combineCleaned(toRemove)=false; %
    [mergedGrains,parentId,combinedTwinBoundary,G.Edges.combineCleaned] = MergeByEdge(G,G.Edges.combineCleaned,grains)    
    %combineMerge comes from CompareBoundaryAndMean
%     G.Edges.combineDiff=G.Edges.combineMerge~=G.Edges.combineCleaned;
    G_clust=rmedge(G,G.Edges.pairs(~G.Edges.combineCleaned,1),...
        G.Edges.pairs(~G.Edges.combineCleaned,2));
    

    
    %Make cases 
    %1) if mean orientation gives grains surrounded by other fragments
    %that are together keep that case 
    %2) Use improvement in convexity measure to break edges (graph paper)
    %3) should we use results to feed back into selection and the final
    %results?
    %4) twin boundary length
    %5) Relaxed Mistol to try getting more innternal grains
    
    %Pull out the global edge labels 
    count=0;
    for i=1:length(G.Edges.pairs)
        if G.Edges.combineCleaned(i)
            count=count+1; 
            G_clust.Edges.GlobalID(count)=i;
        end
    end   
    
    if doplot==true
    %Overlayer clustered grains graph on grains and label edges
        figure;
        plot(grains,grains.meanOrientation,'Micronbar','off','silent');
        hold on 
        plot(mergedGrains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-',...
        'displayName','merged grains')
        plot(combinedTwinBoundary,'linecolor','w','linewidth',2,'displayName','merging boundary');
        p=plot(G_clust,'XData',G_clust.Nodes.centroids(:,1),...
            'YData',G_clust.Nodes.centroids(:,2),'displayName','graph');
        hold off
        p.Marker='s';p.NodeColor='k';p.MarkerSize=3;p.EdgeColor='k';
        if dolabel
            labeledge(p,G_clust.Edges.pairs(:,1),...
                G_clust.Edges.pairs(:,2),G_clust.Edges.GlobalID);
        end
    end
end

