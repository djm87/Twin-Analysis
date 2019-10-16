function [G_clustNew] = FilterTwinBoundaries(G_clust,grains,CS)
%Uses Paris or aspect ratio definition to remove boundaries 
    %Merge all current boundaries by edge
    combine=G_clust.Edges.combineCleaned;
    
    [mergedGrains,parentId,combinedTwinBoundary,combine] = MergeByEdge(G_clust,combine,grains);
    parisAll=parisMerged(mergedGrains,CS);
    
    %Initial variables 
    grainsWithEdges=cell(length(mergedGrains),1);
    maxMergedGrains=1;

    for i=1:length(mergedGrains) 
        grainsWithEdges{i}=grains(parentId == mergedGrains.id(i)).id;
        if length(grainsWithEdges{i})>maxMergedGrains
           maxMergedGrains=maxMergedGrains+1; 
        end
    end   
    
    %Get parent id for edges of all boundaries
    eparentId=parentId(G_clust.Edges.pairs(:,1));
    
    %Initialize variables
    combineKeep=ones(length(G_clust.Edges.pairs),1,'logical');
    %%
    for i=1:maxMergedGrains
        combine=ones(length(G_clust.Edges.pairs),1,'logical');
        for j=1:length(mergedGrains) 
            if length(grainsWithEdges{j})>i
                if ismember(grainsWithEdges{j}(i),G_clust.Edges.pairs) 
                    selected=logical(sum(grainsWithEdges{j}(i)==G_clust.Edges.pairs,2));
                    if sum(selected) == 1
                        combine(selected)=false;
                    end
                end
            end
        end
        combine(~combineKeep)=false;
        [mergedGrainsTest,parentIdTest] = MergeByEdge(G_clust,combine,grains);
        
        %Find parent of id of newly combined edges
        parisTest=parisMerged(mergedGrainsTest,CS);
%         figure;
%         plot(mergedGrainsTest,parisTest,'Micronbar','off','silent');
        
        %Get parent id for edges
        eparentIdTest=parentIdTest(G_clust.Edges.pairs(:,1));

        %for each edge in original and test, compare metrics 
        for j=1:length(eparentIdTest)
           if parisTest(eparentIdTest(j))+0.05*abs(parisAll(eparentId(j)))< parisAll(eparentId(j))
                combineKeep(j)=false;
           end
        end
    end
    G_clustNew=rmedge(G_clust,G_clust.Edges.pairs(~combineKeep,1),...
    G_clust.Edges.pairs(~combineKeep,2));
    combine=G_clustNew.Edges.combineCleaned;
    [mergedGrains] = MergeByEdge(G_clustNew,combine,grains);
        parisFinal=parisMerged(mergedGrains,CS);

        figure;
        plot(mergedGrains,parisFinal,'Micronbar','off','silent');
        hold on 
%         plot(mergedGrains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-',...
%         'displayName','merged grains')
%         plot(combinedTwinBoundary,'linecolor','w','linewidth',2,'displayName','merging boundary');
        p=plot(G_clustNew,'XData',G_clustNew.Nodes.centroids(:,1),...
            'YData',G_clustNew.Nodes.centroids(:,2),'displayName','graph');
         plot(mergedGrains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-',...
        'displayName','merged grains')
        hold off
        p.Marker='s';p.NodeColor='k';p.MarkerSize=3;p.EdgeColor='k';
            labeledge(p,G_clustNew.Edges.pairs(:,1),...
                G_clustNew.Edges.pairs(:,2),G_clustNew.Edges.GlobalID);
        
end

