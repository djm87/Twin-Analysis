function [G_Complete] = filterClusters(G_Complete,grains,useVariantGroup,nVariantGroup)
%filterClusters computes the family matrix
%   Detailed explanation goes here

    for i=1:max(G_Complete.Edges.Group) 
        egroupId = find((i==G_Complete.Edges.Group)==true); %converts logical arrays to indices
        ngroupId = find((i==G_Complete.Nodes.Group)==true);
        nFamily = G_Complete.Nodes.FamilyID(ngroupId);
        nId = G_Complete.Nodes.Id(ngroupId);
        eType = G_Complete.Edges.type(egroupId);
        eVote = G_Complete.Edges.Vote(egroupId,:);
        ePairs = G_Complete.Edges.pairs(egroupId,:);
        eFamily = G_Complete.Edges.FamilyID(egroupId,:);
        eGlobalId = G_Complete.Edges.GlobalID(egroupId);

        [FamilyMatrix,EdgeMatrix,FamilyRelationList,eVotesSummed]=...
             MakeFamilyRelationship(nId,nFamily,eType,eVote,eFamily,ePairs,...
             useVariantGroup,nVariantGroup)
         
        figure(1); 
        plot(grains(nId),...
            G_Complete.Nodes.FamilyID(nId),'Micronbar','off') 
        %{
        figure; 
        plot(grains(nId),...
            grains(nId).meanOrientation,'Micronbar','off')   
                    hold on
        e2keep=(i==G_Complete.Edges.Group)==true;

%             Ggrain=rmedge(G_Complete,G_Complete.Edges.pairs(~e2keep,1),G_Complete.Edges.pairs(~e2keep,2));
        p=plot(G_Complete,'XData',G_Complete.Nodes.centroids(:,1),...
            'YData',G_Complete.Nodes.centroids(:,2),'displayName','graph');
        hold off
        p.Marker='s';p.NodeColor='k';p.MarkerSize=3;p.EdgeColor='k';
        labeledge(p,G_Complete.Edges.pairs(:,1),G_Complete.Edges.pairs(:,2),G_Complete.Edges.GlobalID);
        %}
        G_Complete=G_Complete;
        
        eGlobalId([8,10,11,14,18,19])
    end
end

