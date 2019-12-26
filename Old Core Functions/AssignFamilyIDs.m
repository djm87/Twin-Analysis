function [G_clust,time] = AssignFamilyIDs(G_clust,grains,mergedGrains,seg_angle,twin,doPlot,time)
%Determines the similar orientations that a cluster of fragments share and
%assigns a single family id to those orientations.

    tic

    %Assign each Node a family - assume initially that all nodes are by
    %themselves
    G_clust.Nodes.FamilyID=zeros(length(G_clust.Nodes.Id),1);
    
    %Now handle the graph groups
    FamilyID=cell(max(G_clust.Edges.Group),1);
    typeUnknownLocal=cell(max(G_clust.Edges.Group),1);
    nodeID=cell(max(G_clust.Edges.Group),1);
    Group=G_clust.Nodes.Group;
    typeunknown=G_clust.Nodes.typeUnknown;
    mOri=G_clust.Nodes.meanOrientation;
    parfor i=1:max(G_clust.Nodes.Group)
        %convert logical arrays to indices 
        egroupId= find((i==Group)==true);
        ngroupId= find((i==Group)==true);
        FamilyID{i}=GetFamily(mOri(ngroupId),seg_angle); 
        nodeID{i}=ngroupId;
        typeUnknownLocal{i}=typeunknown(ngroupId);
        %unknow types are ignored during the tree. So that we don't affect the
        %matrix representation of families, make sure unknown types have the
        %highest family numbers 
        if any(typeunknown(ngroupId))
%          %Get unknown and known familes
         unKnownFamily=unique(FamilyID{i}(typeUnknownLocal{i}));
         knownFamily=unique(FamilyID{i}(~typeUnknownLocal{i}));
%          
         for j=1:length(unKnownFamily)
             %Check to see if some of the families have mixed unknown and known
             if any(unKnownFamily(j)==knownFamily)
                 %Assign the relationship of the unknown to known 
                 %This will change how clean Families operates on the
                 %edges
                 typeUnknownLocal{i}(unKnownFamily(j)==FamilyID{i})=true;
             end
         end  
         unKnownFamily=unique(FamilyID{i}(typeUnknownLocal{i}));
         knownFamily=unique(FamilyID{i}(~typeUnknownLocal{i}));
         newFamily=zeros(length(FamilyID{i}),1);
         cnt=1;
         for j=1:length(knownFamily)
             newFamily(knownFamily(j)==FamilyID{i})=cnt;
             cnt=cnt+1;
         end
         for j=1:length(unKnownFamily)
             newFamily(unKnownFamily(j)==FamilyID{i})=cnt;
             cnt=cnt+1;
         end   
         FamilyID{i}=newFamily;
        end
    end
    G_clust.Nodes.typeUnknown(vertcat(nodeID{:}))=vertcat(typeUnknownLocal{:});
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
    
    %Plot the results
    if doPlot==true       
        %Plot edge labeled graph
        labelNodes=false;labelEdges=false;plotG=false;legendOn=false;
        options={'k',5,'s','k',8,'w',2};
        fhandle = plotGraph(grains,mergedGrains,G_clust,...
            G_clust.Nodes.FamilyID,G_clust.Nodes.Id,...
            labelNodes,labelEdges,legendOn,plotG,options);
        
        mtexColorbar;
        mtexTitle('FamilyId')
    end
    
    if ~isfield(time,'AssignFamilyIDs')
        time.AssignFamilyIDs=0;
    end
    time.AssignFamilyIDs=time.AssignFamilyIDs+toc;
end

