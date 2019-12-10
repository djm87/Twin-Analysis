function G = GroupsFromGraph(G,mergedGrains,parentId)
%Given mergedGrains, the function assigns mergedGrains.Id to
%to each edge in the merged grain

    %Initialize
    G.Edges.Group=zeros(length(G.Edges.pairs),1);
    G.Nodes.Group=zeros(length(G.Nodes.Id),1);

    %Asign each pair to a group of grains
    for i=1:length(mergedGrains)
        nId=find(parentId==i);
        G.Nodes.Group(nId)=i;
        for j=1:length(nId)
            G.Edges.Group(nId(j)==G.Edges.pairs(:,1))=i;
        end
    end

end

