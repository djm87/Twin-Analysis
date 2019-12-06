function [G,mergedGrains,mergedGrainsRlx,twinBoundary,twinBoundaryRlx] = MergeByBoundary(G,grains,Mistol,MistolRlx,twin,doplot)
%Applies standard mtex grain merging from boundaries and translates this
%merging to edges where all edges in a merged grain are passed back, not
%just the edges associated with the twin. This is a type of prefilter for
%edge based merging. 

    %Get the ebsd based twin boundary for all twin types in twin
   [twinBoundary] = EBSDTwinBoundary(grains,Mistol,twin);
   [twinBoundaryRlx] = EBSDTwinBoundary(grains,MistolRlx,twin);

    %Combine the twin boundary types together and merge
    combinedTwinBoundary=vertcat(twinBoundary{:});
    [mergedGrains,parentId] = merge(grains,combinedTwinBoundary);
    combinedTwinBoundaryRlx=vertcat(twinBoundaryRlx{:});
    [mergedGrainsRlx,parentIdRlx] = merge(grains,combinedTwinBoundaryRlx);
    
    %Transfer the merged boundaries to edges that were merged
    [combineBoundary] = EdgesFromMergedGrains(mergedGrains,parentId,G.Edges.pairs);
    [combineBoundaryRlx] = EdgesFromMergedGrains(mergedGrainsRlx,parentIdRlx,G.Edges.pairs);
    
    %Create a logical array, where the edge should be kept if true
    G.Edges.combineBoundary=zeros(length(G.Edges.pairs),1,'logical');
    G.Edges.combineBoundary(G.Edges.type==length(twin))=true; %Keep the grain inclusions
    G.Edges.combineBoundary(vertcat(combineBoundary{:}))=true;
    
    G.Edges.combineBoundaryRlx=zeros(length(G.Edges.pairs),1,'logical');
    G.Edges.combineBoundaryRlx(G.Edges.type==length(twin))=true; %Keep the grain inclusions
    G.Edges.combineBoundaryRlx(vertcat(combineBoundaryRlx{:}))=true;

    if doplot
        figure;plot(grains,grains.meanOrientation,'Micronbar','off','silent');    
        hold on 
        plot(mergedGrains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-',...
        'displayName','merged grains')
        colors=hsv(length(twin));
        for i=1:length(twin)-1
            plot(twinBoundary{i},'linecolor',colors(i,:),'linewidth',0.5,'displayName',twin{i}.name);
        end
        
        figure;plot(grains,grains.meanOrientation,'Micronbar','off','silent');    
        hold on 
        plot(mergedGrainsRlx.boundary,'linecolor','k','linewidth',2.5,'linestyle','-',...
        'displayName','merged grains')
        colors=hsv(length(twin));
        for i=1:length(twin)-1
            plot(twinBoundaryRlx{i},'linecolor',colors(i,:),'linewidth',0.5,'displayName',twin{i}.name);
        end        
    end
end

