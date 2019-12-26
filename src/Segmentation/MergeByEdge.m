function [mergedGrains,parentId,combinedEdgeBoundary,GbInd] = MergeByEdge(pairs,GbInd,combine,grains)
%MergeByEdge takes a logical array of edges specifying whether grains
%should be merged. Grain boundaries are extract for those edges and grains
%are grouped accordingly. The out mergedGrains is with respect to the
%original grains. 

    mineral=grains.mineral;
    gB=grains.boundary;
    gB_mineral = gB(mineral,mineral);

    GbIndCombine=GbInd(combine);
    edgeList=find(GbIndCombine==0);
    if ~isempty(edgeList)
        gB_Id=gB_mineral.grainId;
        pairsCombine=pairs(combine,:);  
        %does improve performance Need to rethink how gB are stored/accessed..
        %solution: Only one boundary index is needed for each edge. Simply Save
        %one index so we can merge the edge. light weight in memory and speed.
        for i=1:length(edgeList)
            GbIndCombine(edgeList(i))=find(all(pairsCombine(edgeList(i),:)==gB_Id,2),1,'first');
        end
    end
    GbInd(combine)=GbIndCombine;
     combinedEdgeBoundary=gB_mineral(GbIndCombine);
    [mergedGrains,parentId] = merge(grains,combinedEdgeBoundary);

    
end

