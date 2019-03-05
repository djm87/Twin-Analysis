function [mergedGrains,parentId,combinedTwinBoundary,combine] = MergeByEdge(G,combine,grains)
%MergeByEdge takes a logical array of edges specifying whether grains
%should be merged. Grain boundaries are extract for those edges and grains
%are grouped accordingly. The out mergedGrains is with respect to the
%original grains. 

    mineral=grains.mineral;
    gB=grains.boundary;
    gB_mineral = gB(mineral,mineral);   
    twinBoundary={};
    count=0;
    twinBoundary=cell(sum(combine),1);
    for i=1:length(G.Edges.pairs)   
        if combine(i)
            count=count+1;
            isTwinning=sum(G.Edges.pairs(i,:)==gB_mineral.grainId,2)==2;
            twinBoundary{count} = gB_mineral(isTwinning);
            if size(twinBoundary{count},1)<1
               twinBoundary(count)=[];
               count=count-1; 
               combine(i)=false;
            end
        end
    end

    combinedTwinBoundary=[twinBoundary{:}];
    [mergedGrains,parentId] = merge(grains,combinedTwinBoundary);

end

