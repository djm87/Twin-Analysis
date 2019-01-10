function G_clust = AssignFamilyIDs(G_clust,grains,seg_angle,doplot,dolabel)
%Determines the similar orientations that a cluster of fragments share and
%assigns a single family id to those orientations.

    %remove the grains we don't want to consider in the clustering
    Grains2Keep=unique(G_clust.Edges.pairs);
    G_clust=rmnode(G_clust,G_clust.Nodes.Id(~ismember(G_clust.Nodes.Id,...
        Grains2Keep)));
    %Get the nodes that are connected by edges
    %repeat values mean they are connected
    G_clust.Nodes.Group = conncomp(G_clust)'; 


    %Asign each pair to a group of grains
    G_clust.Edges.Group=zeros(length(G_clust.Edges.pairs),1);
    for i=1:length(G_clust.Edges.pairs)
          G_clust.Edges.Group(i)=G_clust.Nodes.Group(...
              find(G_clust.Edges.pairs(i,1)==G_clust.Nodes.Id));
    end

    %loop over grain clusters to determine  families in grain groups
    G_clust.Nodes.FamilyID=zeros(length(G_clust.Nodes.Id),1);
    for i=1:max(G_clust.Edges.Group)
        %convert logical arrays to indices
        egroupId= find((i==G_clust.Edges.Group)==true); 
        ngroupId= find((i==G_clust.Nodes.Group)==true);

%         [G_clust.Edges.pairs(egroupId,1),G_clust.Edges.pairs(egroupId,2)]
%         G_clust.Nodes.Id(ngroupId)
        
        ori=G_clust.Nodes.meanOrientation(ngroupId);
        G_clust.Nodes.FamilyID(ngroupId)=GetFamily(ori,seg_angle);
        
    end
    
 
    %Determine what family each pair relates
    G_clust.Edges.FamilyID=zeros(length(G_clust.Edges.pairs),2);
    for i=1:max(G_clust.Edges.Group)
        %convert logical arrays to indices
        egroupId= find((i==G_clust.Edges.Group)==true); 
        ngroupId= find((i==G_clust.Nodes.Group)==true);
        nId=G_clust.Nodes.Id(ngroupId);
        fId=G_clust.Nodes.FamilyID(ngroupId);
        for j=1:length(egroupId)
            G_clust.Edges.FamilyID(egroupId(j),1)=...
                unique(fId(G_clust.Edges.pairs(egroupId(j),1)==nId));
            G_clust.Edges.FamilyID(egroupId(j),2)=...
                unique(fId(G_clust.Edges.pairs(egroupId(j),2)==nId));
        end
    end
    if doplot
        figure; 
        plot(grains(G_clust.Nodes.Id),...
            G_clust.Nodes.FamilyID,'Micronbar','off');
        hold on 
        p=plot(G_clust,'XData',G_clust.Nodes.centroids(:,1),...
            'YData',G_clust.Nodes.centroids(:,2));
        hold off
        p.Marker='s';p.NodeColor='k';p.MarkerSize=3;p.EdgeColor='k';
        if dolabel
            labeledge(p,1:length(G_clust.Edges.pairs),...
                1:length(G_clust.Edges.pairs));
        end
    end
end

