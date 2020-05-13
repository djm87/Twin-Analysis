function [combine,mergedGrains,twinBoundary,GBType] = MergeByBoundary(G,grains,mistol,opt)
%Applies standard mtex grain merging from boundaries and translates this
%merging to edges where all edges in a merged grain are passed back, not
%just the edges associated with the twin. This is a type of prefilter for
%edge based merging. 

        

    %Get the ebsd based twin boundary for all twin types in twin
    [twinBoundary] = EBSDTwinBoundary(grains,mistol,opt.nTwin,opt.twin, opt.mergeByGeometry.mergeTP);

    %Construct the merge matrix
    combinedTwinBoundary=vertcat(twinBoundary{:});
    maxId = max(grains.id)+1;
    mergeId=combinedTwinBoundary.grainId;
    M = sparse(mergeId(:,1),mergeId(:,2),1,maxId,maxId);
    [mergedGrains,parentId] = merge(grains,M);

    %Set the grain boundary type and combine flag
    GBType=zeros(size(G.Edges.pairs,1),1,'int8'); 
    combine=zeros(size(G.Edges.pairs,1),1,'logical'); 
    for i=1:opt.nTwin
        GBpairs=unique(twinBoundary{i}.grainId,'rows');
        for j=1:size(GBpairs,1)
            eId=find(all(GBpairs(j,:)==G.Edges.pairs,2) | all(fliplr(GBpairs(j,:))==G.Edges.pairs,2));
            GBType(eId)=i;
            combine(eId)=true;
        end
    end
    
end

