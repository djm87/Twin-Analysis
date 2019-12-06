function [combineBoundary] = EdgesFromMergedGrains(mergedGrains,parentId,pairs)
%Given merged grains find all pairs in each each merged grain and store the
%id of the edges in a cell array of length mergedGrains.
    combineBoundary=cell(length(mergedGrains),1);
    parfor i=1:length(mergedGrains) 
        grainid=find(parentId == mergedGrains(i).id);
        if length(grainid)>1
            combine1=zeros(length(pairs(:,1)),1,'logical');
            combine2=zeros(length(pairs(:,1)),1,'logical');
            
            %for each merged grain, identify all edges that have the grains
            for j=1:length(grainid)
                combine1(find(grainid(j)==pairs(:,1)))=true;
                combine2(find(grainid(j)==pairs(:,2)))=true;
            end
            %Sum the pairs to see which edges are in the merged
            %grains. This creates a combine that is based on boundary
            combineBoundary{i}=find((combine1+combine2)==2);
        end
    end
end