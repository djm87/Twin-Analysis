function [combine,mergedGrains,twinBoundary,GBType] = MergeByBoundary(G,grains,mistol,opt)
%Applies standard mtex grain merging from boundaries and translates this
%merging to edges where all edges in a merged grain are passed back, not
%just the edges associated with the twin. This is a type of prefilter for
%edge based merging. 

        

    %Get the ebsd based twin boundary for all twin types in twin
    [twinBoundary] = EBSDTwinBoundary(grains,mistol,opt.nMori,opt.mori, opt.gclust.mergeTP);

    %Construct the merge matrix
    combinedTwinBoundary=vertcat(twinBoundary{:});
    [mergedGrains,parentId] = merge(grains,combinedTwinBoundary);
    
    %Set the grain boundary type and combine flag
    GBType=zeros(size(G.Edges.pairs,1),1,'int8'); 
    combine=zeros(size(G.Edges.pairs,1),1,'logical'); 
    for i=1:opt.nMori
        GBpairs=unique(twinBoundary{i}.grainId,'rows');
        for j=1:size(GBpairs,1)
            eId=find(all(GBpairs(j,:)==G.Edges.pairs,2) | all(fliplr(GBpairs(j,:))==G.Edges.pairs,2));
            GBType(eId)=i;
            combine(eId)=true;
        end
    end
    
end

