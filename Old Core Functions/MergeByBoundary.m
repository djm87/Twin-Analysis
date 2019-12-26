function [combineBoundary,mergedGrains,twinBoundary,time] = MergeByBoundary(G,grains,Mistol,twin,doplot,time)
%Applies standard mtex grain merging from boundaries and translates this
%merging to edges where all edges in a merged grain are passed back, not
%just the edges associated with the twin. This is a type of prefilter for
%edge based merging. 

    tic

    %Get the ebsd based twin boundary for all twin types in twin
   [twinBoundary] = EBSDTwinBoundary(grains,Mistol,twin);

    %Combine the twin boundary types together and merge
    combinedTwinBoundary=vertcat(twinBoundary{:});
    [mergedGrains,parentId] = merge(grains,combinedTwinBoundary);

    
    %Transfer the merged boundaries to edges that were merged
    [combineBoundaryTmp] = EdgesFromMergedGrains(mergedGrains,parentId,G.Edges.pairs);
    
    %Create a logical array, where the edge should be kept if true
    combineBoundary=zeros(length(G.Edges.pairs),1,'logical');
    combineBoundary(G.Edges.type==length(twin))=true; %Keep the grain inclusions
    combineBoundary(vertcat(combineBoundaryTmp{:}))=true;

    if doplot
        figure;plot(grains,grains.meanOrientation,'Micronbar','off','silent');    
        hold on 
        plot(mergedGrains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-',...
        'displayName','merged grains')
        colors=hsv(length(twin));
        for i=1:length(twin)-1
            plot(twinBoundary{i},'linecolor',colors(i,:),'linewidth',0.5,'displayName',twin{i}.name);
        end
        hold off
        mtexTitle('Combined with boundary')
        legend off

    end
    if ~isfield(time,'MergeByBoundary')
        time.MergeByBoundary=0;
    end
    time.MergeByBoundary=time.MergeByBoundary+toc;
end

