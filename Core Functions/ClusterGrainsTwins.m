function G_clust = ClusterGrainsTwins(G,grains,Mistol,meanMistolRelaxed,twin,...
    removeList,addList,doplot,dolabel,doTwinGrainOnly)
%Clusters grains by removing edges of grains that should not be combined
    %Returns G.Edges.combineBoundary and G.Edges.groupMerge (currently not
    %returned)
    
    [G,mergedGrains,twinBoundary]= MergeByBoundary(G,grains,Mistol,twin);
    
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
    G.Edges.combineCleaned(addList)=true;
    G.Edges.combineCleaned(toRemove)=false; %
    [mergedGrains,parentId,combinedTwinBoundary,G.Edges.combineCleaned] = MergeByEdge(G,G.Edges.combineCleaned,grains);    
    %combineMerge comes from CompareBoundaryAndMean
%     G.Edges.combineDiff=G.Edges.combineMerge~=G.Edges.combineCleaned;
    G_clust=rmedge(G,G.Edges.pairs(~G.Edges.combineCleaned,1),...
        G.Edges.pairs(~G.Edges.combineCleaned,2));
    G_clust_removed=rmedge(G,G.Edges.pairs(G.Edges.combineCleaned,1),...
        G.Edges.pairs(G.Edges.combineCleaned,2));   
    
    %Pull out the global edge labels 
    count=0;
    count_removed=0;
    for i=1:length(G.Edges.pairs)
        if G.Edges.combineCleaned(i)
            count=count+1; 
            G_clust.Edges.GlobalID(count)=i;
        else
            count_removed=count_removed+1; 
            G_clust_removed.Edges.GlobalID(count_removed)=i;
        end
    end 
    

    
    if doplot==true
        %Remove unused nodes so that matlab plot performance is better
        toremove=ones(length(G_clust.Nodes.Id),1,'logical');
        toremove(unique(G_clust.Edges.pairs))=false;
        G_clust_small=G_clust;
        G_clust_small=rmnode(G_clust,find(toremove));
        %Overlayer clustered grains graph on grains and label edges
        figure; 
        if doTwinGrainOnly
            plot(grains(G_clust_small.Nodes.Id),...
                G_clust_small.Nodes.meanOrientation,'Micronbar','off','silent');
        else
            plot(grains(G_clust.Nodes.Id),...
                G_clust.Nodes.meanOrientation,'Micronbar','off','silent');
        end
        hold on 
%         plot(mergedGrains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-',...
%         'displayName','merged grains')
        colors=hsv(length(twin));
        for i=1:length(twin)
            plot(twinBoundary{i},'linecolor',colors(i,:),'linewidth',0.5,'displayName',twin{i}.name);
        end
%         plot(combinedTwinBoundary,'linecolor','w','linewidth',2,'displayName','merging boundary');
        p=plot(G_clust_small,'XData',G_clust_small.Nodes.centroids(:,1),...
            'YData',G_clust_small.Nodes.centroids(:,2),'displayName','graph');
         plot(mergedGrains.boundary,'linecolor','k','linewidth',0.5,'linestyle','-',...
        'displayName','merged grains')
        hold off
        p.EdgeColor='k';p.MarkerSize=5;p.Marker='s';p.NodeColor='k';
        if dolabel
            pairs1=G_clust_small.Edges.pairs(:,1);
            pairs2=G_clust_small.Edges.pairs(:,2);
            for i=1:length(G_clust_small.Nodes.Id)
                pairs1(pairs1==G_clust_small.Nodes.Id(i))=i;
                pairs2(pairs2==G_clust_small.Nodes.Id(i))=i;
            end
            labeledge(p,pairs1,...
                pairs2,G_clust_small.Edges.GlobalID);
%         p.EdgeFontSize=1.5;
        legend off
%         print('Edges','-dpdf','-r600');
        end 
    end
end

