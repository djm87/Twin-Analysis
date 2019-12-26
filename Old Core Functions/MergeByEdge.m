function [mergedGrains,parentId,combinedTwinBoundary,combine] = MergeByEdge(pairs,combine,grains)
%MergeByEdge takes a logical array of edges specifying whether grains
%should be merged. Grain boundaries are extract for those edges and grains
%are grouped accordingly. The out mergedGrains is with respect to the
%original grains. 
    mineral=grains.mineral;
    gB=grains.boundary;
    gB_mineral = gB(mineral,mineral);
    gB_Id=gB_mineral.grainId;
    pairsCombine=pairs(combine,:);

    isTwinBoundary=zeros(length(gB_Id),1,'logical');
    for i=1:length(pairsCombine)   
        isTwinBoundary(all(pairsCombine(i,:)==gB_Id,2)) = true;
    end

    combinedTwinBoundary=gB_mineral(isTwinBoundary);
    [mergedGrains,parentId] = merge(grains,combinedTwinBoundary);

end

