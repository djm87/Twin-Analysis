function [G_clust,G] = AssignFamilyIDs(G_clust,G,grains,mGrains,opt)
%Determines the similar orientations that a cluster of fragments share and
%assigns a single family id to those orientations.
    groupList=unique(G_clust.Nodes.Group(G_clust.Nodes.isNewGroup));
    if isempty(groupList)
       return 
    end
    
    %Now handle the graph groups
    FamilyID=cell(length(groupList),1);
    typeUnknownLocal=cell(length(groupList),1);
    nodeID=cell(length(groupList),1);
    groups=G_clust.Nodes.Group;
    typeunknown=G_clust.Nodes.typeUnknown;
    mOri=grains.meanOrientation;
    for i=1:length(groupList)
        group=groupList(i);
        %convert logical arrays to indices 
        egroupId= find((group==groups));
        ngroupId= find((group==groups));

        FamilyID{i}=GetFamily(mOri(ngroupId),opt.grain_recon.seg_angle_grouped); 
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
    G.Nodes.FamilyID(:)=0; %reset
    G.Nodes.FamilyID(G_clust.Nodes.Id)=G_clust.Nodes.FamilyID;

    G_clust.Edges.FamilyID(:,1)=G_clust.Nodes.FamilyID(G_clust.Edges.pairs(:,1));
    G_clust.Edges.FamilyID(:,2)=G_clust.Nodes.FamilyID(G_clust.Edges.pairs(:,2));
    G.Edges.FamilyID(:,:)=0; %reset
    G.Edges.FamilyID(G.Edges.combineCleaned,:)=G_clust.Edges.FamilyID;

    %Plot the results
    if opt.plot.do    
        %Plot edge labeled graph
        labelNodes=false;labelEdges=opt.plot.labelEdges;plotG=false;legendOn=opt.plot.legendOn;
        options={'k',5,'s','k',8,'w',2};
        fhandle = plotGraph(grains,mGrains,G_clust,...
            G_clust.Nodes.FamilyID,G_clust.Nodes.Id,...
            labelNodes,labelEdges,legendOn,plotG,options);
        mtexColorbar;
        mtexTitle('FamilyId')
    end
    
end

