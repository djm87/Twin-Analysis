function [twinCount,G_Complete] = CountTwins(G_Complete)
%CountTwins produces twin count statistics 
    twinCount = zeros(max(G_Complete.Edges.Group) ,1);
    G_Complete.Nodes.twinCount=zeros(length(G_Complete.Nodes.FamilyID) ,1);
    for i=1:max(G_Complete.Edges.Group) 
        egroupId = find((i==G_Complete.Edges.Group)==true); %converts logical arrays to indices
        ngroupId = find((i==G_Complete.Nodes.Group)==true);
        nMergeTwin = G_Complete.Nodes.MergeTwin(ngroupId);
        nType = G_Complete.Nodes.Type(ngroupId); 
        nId = G_Complete.Nodes.Id(ngroupId);

        %Find the unique twins to merge
        twins = unique(nMergeTwin);
        
        %If nMergeTwin>0 only count the merged twins once + only count twin
        %type fragments.
        twinCount(i)=length(twins)-1+sum(nMergeTwin==0 & nType>0);
        G_Complete.Nodes.twinCount(nId)=twinCount(i);
    end
end

