function [nGroup,eGroup] = GroupsFromGraph(pairs,nGrains,nMergedGrains,parentId)
%Given mergedGrains, the function assigns mergedGrains.Id to
%to each edge in the merged grain

    %Initialize
    eGroup=zeros(size(pairs,1),1);
    nGroup=zeros(nGrains,1);

    %Asign each pair to a group of grains
    for i=1:nMergedGrains
        nId=find(parentId==i);
        nGroup(nId)=i;
        [~,loc]=intersect(pairs(:,1),nId,'rows');
        eGroup(loc)=i;
    end

end

