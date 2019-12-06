function [G_clust,mergedGrains] = ClusterGrainsTwins(G,grains,meanMistol,meanMistolRelaxed,minNEdgeMistol,twin,...
    eRemoveList,eAddList,nRemoveList,nAddList,doplot,dolabel,doTwinGrainOnly)
    %Clusters grains by removing edges of grains that should not be combined
    %Returns G.Edges.combineBoundary 
    %Initialize globalId
    G.Edges.GlobalID=[1:length(G.Edges.pairs(:,1))]';
    
    %Compute misorientation between all pairs
    mori=inv(grains(G.Edges.pairs(:,1)).meanOrientation).*...
        grains(G.Edges.pairs(:,2)).meanOrientation; 

    %Get the type of misorientation
    [~,type] = TestTwinRelationship(mori,meanMistol,twin,G.Edges.type);
    [~,typeRlx] = TestTwinRelationship(mori,meanMistolRelaxed,twin,G.Edges.type);

    G.Edges.type=typeRlx;
    twinEdges=G.Edges.type > 0;
    
    %Any edges that aren't a twin should be set to be removed
    %Filter edge by mistol based with meanMistolRelaxed
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
    
    %Add internal 
    twinEdgeNodes=unique(G.Edges.pairs([twinEdges,twinEdges]));
    G.Nodes.typeUnknown(twinEdgeNodes)=false;
    
    %Try removing nodes
    for i=1:length(nRemoveList)
        [~,epairs]=grains(nRemoveList(i)).neighbors;
        for j=1:size(epairs,1)
            eId=find(all(epairs(j,:)==G.Edges.pairs,2) | all(fliplr(epairs(j,:))==G.Edges.pairs,2));
            isTwin=G.Edges.type(eId)~=0 & G.Edges.type(eId)~=length(twin);
            G.Edges.combineCleaned(eId(isTwin))=false;
        end
    end
    
    %Remove edges
    G.Edges.combineCleaned(eRemoveList)=false; %

    %Add edges
    typeAdd=typeRlx(eAddList);
    typeAdd(typeAdd==0)=length(twin); %Unknown Type
    G.Edges.type(eAddList)=typeAdd;
    G.Edges.combineCleaned(eAddList)=true;

    %Update typeUnknown list with current edge types (can change from
    %initial and different mergedGrain configurations

    
    [mergedGrains,parentId,combinedTwinBoundary,G.Edges.combineCleaned] = MergeByEdge(G.Edges.pairs,G.Edges.combineCleaned,grains,typeRlx,minNEdgeMistol);    
    
    %DS: There 
%     %label internal grains as type length(twin) 
%     isInside = checkInside(mergedGrains, mergedGrains);
%     [GrainIdInclusion,GrainIdWithInclusion] = find(isInside);
%     figure;plot(mergedGrains(GrainIdWithInclusion),mergedGrains(GrainIdWithInclusion).area)
%     for i=1:length(GrainIdInclusion)
%         grainsInclusion=grains(parentId==mergedGrains(GrainIdInclusion(i)).id);
%         G.Nodes.typeUnknown(grainsInclusion.id)=true;
% 
%         %Apply edges to all neighbors (these factor into the family tree)        
%         [~,epairs]=grainsInclusion.neighbors;
%         for j=1:size(epairs,1)
%             eId=find(all(epairs(j,:)==G.Edges.pairs,2) | all(fliplr(epairs(j,:))==G.Edges.pairs,2));
%             isTwin=G.Edges.type(eId)~=0 & G.Edges.type(eId)~=length(twin);
%             if any(isTwin)
%                G.Edges.combineCleaned(eId(isTwin))=true;
%             else
%                 G.Edges.type(eId)=length(twin);
%                 G.Edges.combineCleaned(eId)=true;
%             end
%         end
%     end

    %Check if internal grains have a type from relaxed tol
    

    
    
    %Remerge Grains for visualization
%     [mergedGrains,parentId,combinedTwinBoundary,G.Edges.combineCleaned] = MergeByEdge(G.Edges.pairs,G.Edges.combineCleaned,grains);    
    
    %combineMerge comes from CompareBoundaryAndMean
%     G.Edges.combineDiff=G.Edges.combineMerge~=G.Edges.combineCleaned;
%     mergedInteralGrains=mergedGrains.inclusionId
%     GrainWithInternal=find(mergedInteralGrains>0);
%     figure;plot(grains(GrainWithInternal),grains(GrainWithInternal).meanOrientation)
%     hold on;
%     plot(mergedGrains.boundary)
%     hold off
%     for i=1:length(mergedGrains)
%         mergedGrainId=find(mergedInteralGrains(i,:));
%         if ~isempty(mergedGrainId)
%             disp('hello')
%             for j=mergedGrainId
%                 grains.neighbors    
%             end
%         end
%     end
    
    G_clust=rmedge(G,G.Edges.pairs(~G.Edges.combineCleaned,1),...
        G.Edges.pairs(~G.Edges.combineCleaned,2));
%     G_clust_removed=rmedge(G,G.Edges.pairs(G.Edges.combineCleaned,1),...
%         G.Edges.pairs(G.Edges.combineCleaned,2));   
    
%     [uparentId i j] = unique(parentId,'first');
%     indexToDupes = find(not(ismember(1:numel(parentId),i)))
    %Pull out the global edge labels 
%     count=0;
%     count_removed=0;
%     for i=1:length(G.Edges.pairs)
%         if G.Edges.combineCleaned(i)
%             count=count+1; 
%             G_clust.Edges.GlobalID(count)=i; %Just assign this before removing edges?
%         else
%             count_removed=count_removed+1; 
%             G_clust_removed.Edges.GlobalID(count_removed)=i;
%         end
%     end 
    

    
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

        p=plot(G_clust_small,'XData',G_clust_small.Nodes.centroids(:,1),...
            'YData',G_clust_small.Nodes.centroids(:,2),'displayName','graph');
         plot(mergedGrains.boundary,'linecolor','k','linewidth',2,'linestyle','-',...
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
        
        figure; 
        
        if doTwinGrainOnly
            plot(grains(G_clust_small.Nodes.Id),...
                G_clust_small.Nodes.meanOrientation,'Micronbar','off','silent');
        else
            plot(grains(G_clust.Nodes.Id),...
                G_clust.Nodes.meanOrientation,'Micronbar','off','silent');
        end
        
        hold on
         plot(mergedGrains.boundary,'linecolor','k','linewidth',1.5,'linestyle','-',...
            'displayName','merged grains')
        text(grains(G_clust.Nodes.Id),int2str(G_clust.Nodes.Id))
        hold off
        
    end
end

