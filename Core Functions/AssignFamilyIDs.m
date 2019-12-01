function G_clust = AssignFamilyIDs(G_clust,grains,seg_angle,doplot,dolabel,doTwinGrainOnly)
%Determines the similar orientations that a cluster of fragments share and
%assigns a single family id to those orientations.

    %remove the grains we don't want to consider in the clustering
    Grains2Keep=unique(G_clust.Edges.pairs);
%     G_clust=rmnode(G_clust,G_clust.Nodes.Id(~ismember(G_clust.Nodes.Id,...
%         Grains2Keep)));
    %Get the nodes that are connected by edges
    %repeat values mean they are connected
    G_clust.Nodes.Group = conncomp(G_clust)'; 
    
    %Asign each pair to a group of grains
    G_clust.Edges.Group=zeros(length(G_clust.Edges.pairs),1);
    for i=1:length(G_clust.Edges.pairs)
          G_clust.Edges.Group(i)=G_clust.Nodes.Group(...
              find(G_clust.Edges.pairs(i,1)==G_clust.Nodes.Id));
    end
    
    %Assign each Node a family - assume initially that all nodes are by
    %themselves
    G_clust.Nodes.FamilyID=zeros(length(G_clust.Nodes.Id),1);
    
    %Now handle the graph groups
    FamilyID=cell(max(G_clust.Edges.Group),1);
    nodeID=cell(max(G_clust.Edges.Group),1);
    Group=G_clust.Nodes.Group;
    mOri=G_clust.Nodes.meanOrientation;
    parfor i=1:max(G_clust.Nodes.Group)
        %convert logical arrays to indices 
        ngroupId= find((i==Group)==true);
%         FamilyID(ngroupId)=GetFamily(mOri(ngroupId),seg_angle);   
        FamilyID{i}=GetFamily(mOri(ngroupId),seg_angle); 
        nodeID{i}=ngroupId; 
    end
    G_clust.Nodes.FamilyID(vertcat(nodeID{:}))=vertcat(FamilyID{:});
 
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
        %Remove unused nodes so that matlab plot performance is better
        toremove=ones(length(G_clust.Nodes.Id),1,'logical');
        toremove(unique(G_clust.Edges.pairs))=false;
        G_clust_small=rmnode(G_clust,find(toremove));
        
        figure; 
        if doTwinGrainOnly
            
            plot(grains(G_clust_small.Nodes.Id),...
                G_clust_small.Nodes.FamilyID,'Micronbar','off');
        else
            plot(grains(G_clust.Nodes.Id),...
                G_clust.Nodes.FamilyID,'Micronbar','off');
        end
        hold on 
        p=plot(G_clust_small,'XData',G_clust_small.Nodes.centroids(:,1),...
            'YData',G_clust_small.Nodes.centroids(:,2));
        hold off

        p.Marker='s';p.NodeColor='k';p.MarkerSize=3;p.EdgeColor='k';
        if dolabel
            pairs1=G_clust_small.Edges.pairs(:,1);
            pairs2=G_clust_small.Edges.pairs(:,2);
            for i=1:length(G_clust_small.Nodes.Id)
                pairs1(pairs1==G_clust_small.Nodes.Id(i))=i;
                pairs2(pairs2==G_clust_small.Nodes.Id(i))=i;
            end
            labeledge(p,pairs1,...
                pairs2,G_clust_small.Edges.GlobalID);
        end
    end
end

