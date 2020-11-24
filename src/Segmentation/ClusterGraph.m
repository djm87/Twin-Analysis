function [G_clust,G_frag,G_family,mGrains] = ClusterGraph(G_frag,G_family,grains,opt,recompute)
    %ClusterGraph performs the grouping based on edge definitions and user input
    %files. In addition, the families in each cluster are computed and 
    %stored for later use.
        
    %Using the twin Edge list initialized in the G_frag, user defined list, and
    %geometry merging defined a list of edges to keep i.e. G_frag.Edges.combineCleaned
    %In addition we need to make sure the type of relationship i.e. twin,
    %unknown, etc... are appropriately defined
    %======================================================================
    G_frag=getEdgeList(G_frag,opt);
    
    %Remove the edges to compute G_clust and compute the merged grains
    %======================================================================
    G_clust=rmedge(G_frag,G_frag.Edges.pairs(~G_frag.Edges.combineCleaned,1),...
        G_frag.Edges.pairs(~G_frag.Edges.combineCleaned,2));

    [mGrains,parentId] = merge(grains,G_clust.Edges.pairs);    
    
    %Update which cluster has changed compared to the previously computed
    %G_clust that is stored in G_frag and reset properties if needed
    %======================================================================
    %Assign group Ids using mGrains parentId
    G_clust.Edges.Group=parentId(G_clust.Edges.pairs(:,1));
    G_clust.Nodes.Group=parentId;
    G_frag.Nodes.Group=parentId;
    
    %Transfer still valid Family and family vote info to G_clust.. build
    %list of groups that need to be computed
    [G_clust,G_frag,G_family] = GetGroupsToCompute(G_clust,G_frag,G_family,recompute);
    
    %Update G_frag with groups in G_clust
    G_frag.Nodes.Group(G_clust.Nodes.Id)=G_clust.Nodes.Group;
    G_frag.Edges.Group(:)=0;
    G_frag.Edges.Group(G_frag.Edges.combineCleaned)=G_clust.Edges.Group;
    
    %Update properties for grains that need computing
    %======================================================================
    %Update unknown nodes
    nId_edges=unique(G_clust.Edges.pairs([G_clust.Edges.type~=opt.moriUnknown,G_clust.Edges.type~=opt.moriUnknown]));
    G_clust.Nodes.typeUnknown(:)=true; 
    G_clust.Nodes.typeUnknown(nId_edges)=false;

    % Identify Families in grain clusters 
    [G_clust,G_frag]= AssignFamilyIDs(G_clust,G_frag,grains,mGrains,opt);
    
    % Enforce clusters with single Family to have type = 0
    % This is enforced because types are saved onced the Family has been computed and comes back
    % clean
    [G_clust,G_frag] = SetTypeForUntwinnedCluster(G_clust,G_frag);
    
    % Reset the isNewGroup flag to prevent
    G_clust.Nodes.isNewGroup(:)=false; 
    G_frag.Nodes.isNewGroup(:)=false;
end

function [G_frag]=getEdgeList(G_frag,opt)

    % Test if all the neccessary files exist and generate if they don't
    flist={'eAddList.txt','eRemoveList.txt','nRemoveList.txt','nAddList.txt'};
    generateInputs(flist); %move to cluster
    
    % Load content of files
    eAddList=load(flist{1}); %Adds edge and tries to give it a twin relationship, otherwise type twin unknown
    eAddList(find(any(eAddList>G_frag.numedges,2)),:)=[];
    eRemoveList=load(flist{2}); %Remove the edge 
    eRemoveList(eRemoveList>G_frag.numedges)=[];
    nRemoveList=load(flist{3}); %Removes all edges associated with node
    nAddList=load(flist{4}); %Tries to add all edges that have twin relationship
        
    %As a starting point, use tight boundary misorientation. Both
    %combine and GBType are not changed after initialization
    G_frag.Edges.combineCleaned = G_frag.Edges.combine;
    G_frag.Edges.type=G_frag.Edges.typeIni;
    
    %To simplify edge id extraction from a node pair
    pairs=G_frag.Edges.pairs; 
    maxEdgeId=size(pairs,1);
    AllPairs=vertcat(pairs,fliplr(pairs));
    AllPairId=double(vertcat(G_frag.Edges.GlobalID,G_frag.Edges.GlobalID));
    eEdgeId=sparse(AllPairs(:,1),AllPairs(:,2),AllPairId,maxEdgeId,maxEdgeId);
    
    %Try adding all edges attached to a user specified node using a relaxed tolerance
    eToAdd=nonzeros(eEdgeId(nAddList,:));
    eToAdd=eToAdd(G_frag.Edges.combineRlx(eToAdd));
    G_frag.Edges.combineCleaned(eToAdd)=true;
    G_frag.Edges.type(eToAdd)=G_frag.Edges.GBTypeRlx(eToAdd);
    
    %Add user specified edge
    eToAdd=zeros(size(eAddList,1),1);
    for i=1:size(eAddList,1)
        tmp=nonzeros(eEdgeId(eAddList(i,1),eAddList(i,2)));
        if ~isempty(tmp)
            eToAdd(i)=tmp;
        end
    end
    eToAdd(eToAdd==0)=[];
    eTypeAdd=G_frag.Edges.meanType(eToAdd);
    eTypeAdd(eTypeAdd==0)=opt.moriUnknown; %unknown type
    G_frag.Edges.type(eToAdd)=eTypeAdd;
    G_frag.Edges.combineCleaned(eToAdd)=true;
    
    %Update unknown types in nodes
    nId_edges=unique(G_frag.Edges.pairs([G_frag.Edges.type~=opt.moriUnknown,G_frag.Edges.type~=opt.moriUnknown]));
    G_frag.Nodes.typeUnknown(:)=true;
    G_frag.Nodes.typeUnknown(nId_edges)=false;
    
    %Removing all edges attached to a user specified node
    toRemove=nonzeros(eEdgeId(nRemoveList,:));
    G_frag.Edges.combineCleaned(toRemove)=false;
    
    %remove user specified edge
    eRemoveList(eRemoveList>G_frag.numedges)=[];
    G_frag.Edges.combineCleaned(eRemoveList)=false;
    
    %Add back the geometry based merging when there is some geometry to
    %test it on
    if opt.gclust.mergeCluster
        G_frag = MergeByGeometry(G_frag,eEdgeId,grains,opt)
    end
    
    %Update unknown types in nodes
    nId_edges=unique(G_frag.Edges.pairs([G_frag.Edges.type~=opt.moriUnknown,G_frag.Edges.type~=opt.moriUnknown]));
    G_frag.Nodes.typeUnknown(:)=true;
    G_frag.Nodes.typeUnknown(nId_edges)=false;
    
    %Removing all edges attached to a user specified node
    toRemove=nonzeros(eEdgeId(nRemoveList,:));
    G_frag.Edges.combineCleaned(toRemove)=false;
    
    %remove user specified edge
    G_frag.Edges.combineCleaned(eRemoveList)=false;
    
    %Ensure no zero types (This should always be true unless there is an error in logic).
    assert(all(G_frag.Edges.type(G_frag.Edges.combineCleaned)~=0),'There is an error in type definition in clusterGraph')
end

function [G_clust,G_frag,G_family] = GetGroupsToCompute(G_clust,G_frag,G_family_old,recompute)
    %Get group ids to process
    %New clusters should have properties recomputed and old clusters
    %shouldn't. To this end, when a cluster is new relevant properties are
    %reset, a new group flag is set so that Families are computed in 
    %ClusterGraph,  a flag is set for FamilyGraph computation, and
    %FamilyTree computation.
    
    %Get the number of unique groups in G_clust
    nGroup=max(G_clust.Nodes.Group);
    
    %Initialize the G_Family array of graphs
    
    G_family=cell(nGroup,1);
    if isempty(G_family_old) G_family_old=G_family;end
        
    if recompute
        %Recompute everything
        G_clust.Nodes.isNewGroup(:)=true;
        G_clust.Nodes.computeGfam(:)=true;
        G_clust.Nodes.computeTree(:)=true;
    else
        for i=1:nGroup
            %For a given G_clust cluster get one grain
            eID=G_clust.Edges.GlobalID(i==G_clust.Edges.Group);
            nID = find((i==G_clust.Nodes.Group));

            %Find the stored group in G_frag that contains that grain
            oldGroupID=G_frag.Nodes.Group(nID(1));
            G_nID = find(oldGroupID==G_frag.Nodes.Group);
            G_eId = find(oldGroupID==G_frag.Edges.Group);

            %If the store node ids and edge ids match the node ids of G_clust, 
            %then the cluster has not changed
            if length(G_nID) ~= length(nID) || length(G_eId) ~= length(eID) 
    %             groupCompute(i)=true;
                G_clust.Nodes.isNewGroup(nID)=true;
                G_clust.Nodes.computeTree(nID)=true;
            elseif ~all(nID==G_nID) || ~all(eID==G_eId) %Because sorted
    %             groupCompute(i)=true; 
                G_clust.Nodes.isNewGroup(nID)=true;
                G_clust.Nodes.computeTree(nID)=true;
            elseif ~isempty(G_family_old{oldGroupID})
                %Repopulate G_Family to correspond with the new groups
                G_family(i)=G_family_old{oldGroupID};
            end
        end
    end
%     groupList=find(groupCompute);
    %Store in persistent graph
    G_frag.Nodes.isNewGroup(G_clust.Nodes.Id)=G_clust.Nodes.isNewGroup;
    G_frag.Nodes.computeTree(G_clust.Nodes.Id)=G_clust.Nodes.computeTree;

    %Reset the Family information
    G_clust.Nodes.Generation(G_frag.Nodes.isNewGroup)=-1;
    G_frag.Nodes.Generation=G_clust.Nodes.Generation;
    G_clust.Nodes.type(G_frag.Nodes.isNewGroup)=0;
    G_frag.Nodes.type=G_clust.Nodes.type;
    G_clust.Nodes.ParentFamilyId(G_frag.Nodes.isNewGroup)=0;
    G_frag.Nodes.ParentFamilyId =G_clust.Nodes.ParentFamilyId;
    G_clust.Nodes.SchmidVariant(G_frag.Nodes.isNewGroup)=0;
    G_frag.Nodes.SchmidVariant=G_clust.Nodes.SchmidVariant;
    G_clust.Nodes.K1NAng(G_frag.Nodes.isNewGroup)=0;
    G_frag.Nodes.K1NAng=G_clust.Nodes.K1NAng;
    
end

function [G_frag] = MergeByGeometry(G_frag,eEdgeId,grains,opt)
       
    %Initialize variables
    ori=grains.meanOrientation;
%     isGNeighbor=zeros(length(grains),length(grains),'logical');
%     ind = sub2ind(size(isGNeighbor), ePairs(:,1), ePairs(:,2));
%     isGNeighbor(ind)=true;

%     pairs=G_frag.Edges.pairs; 
%     maxEdgeId=size(pairs,1);
%     AllPairs=vertcat(pairs,fliplr(pairs));
%     AllPairId=double(vertcat(G_frag.Edges.GlobalID,G_frag.Edges.GlobalID));
%     eEdgeId=sparse(AllPairs(:,1),AllPairs(:,2),AllPairId,maxEdgeId,maxEdgeId);
    
    ePairs=G_frag.Edges.pairs;
    cnt=1;
    while cnt<=opt.mergeByGeometry.itter    
        %Merge grains based on edge combineCleaned list
        [mGrains,parentId,mTwinGb,G_frag.Edges.GbInd] = MergeByEdge(G_frag.Edges.pairs,G_frag.Edges.GbInd,G_frag.Edges.combineCleaned,grains);    
        
        %Rebuild the family mean orientations 
        mFamilyCenters=cell(length(mGrains),1);
        for i=1:length(mGrains)
            [~,mFamilyCenters{i}]=calcCluster(ori(parentId==i),'maxAngle',opt.mergeByGeometry.FamilyMisTol,'method','hierarchical');
%             misFamily=unique(angle_outer(mFamilyCenters{3214},mFamilyCenters{3214})./degree)
%             min(misFamily(misFamily>opt.grain_recon.FamilyMisTol))/2;
        end
        
        %Plot full graph
%         labelNodes=true;labelEdges=false;plotG=false;legendOn=false;
%         fhandle = plotGraph(grains,mGrains,G_frag,...
%          grains.meanOrientation,1:length(grains),...
%          labelNodes,labelEdges,legendOn,plotG,[]);
        
        %Get the ratios of boundary shared by merged grains
        SBR = SharedBoundaryRatio(mGrains);

        %Create merged graph to handle cluster merging
        [~,mPairs]=neighbors(mGrains);
        s=mPairs(:,1);
        t=mPairs(:,2);
        G_geo=graph(s,t);
        G_geo.Edges.GlobalID=[1:G_geo.numedges]';
        
        %To simplify edge id extraction from a node pair
        mPairs=G_geo.Edges.EndNodes; 
        maxEdgeId=size(mPairs,1);
        AllmPairs=vertcat(mPairs,fliplr(mPairs));
        AllmPairId=double(vertcat(G_geo.Edges.GlobalID,G_geo.Edges.GlobalID));
        meEdgeId=sparse(AllmPairs(:,1),AllmPairs(:,2),AllmPairId,maxEdgeId,maxEdgeId);
        
        %Get small grains for merging
        G_geo.Nodes.small=mGrains.grainSize< opt.mergeByGeometry.SGT; 
        
        %Get high aspect ratio grains for merging 
        G_geo.Nodes.largeAR=mGrains.aspectRatio > opt.mergeByGeometry.ART &...
            mGrains.grainSize < opt.mergeByGeometry.MARGST ; 
        
%         figure;plot(grains([nId_clust;nId]),grains([nId_clust;nId]).meanOrientation)
% %         figure;plot(grains([nId_clust]),grains([nId_clust]).meanOrientation)
%         figure;plot(mGrains([510,9]),mGrains([510,9]).area)

        %To prevent
        meGID=zeros(G_geo.numedges,1);
        mnType=zeros(G_geo.numnodes,1);
        mnSBR=zeros(G_geo.numnodes,1,'logical');
        
        combineMergeByGeometry=zeros(size(ePairs,1),1,'logical');
        GBTypeRlx=G_frag.Edges.GBTypeRlx;
        pairs=G_frag.Edges.pairs;
        for i=1:length(mGrains)
            %Try merging based on boundary ratios. i.e. if cluster/fragment
            %is mostly internal to another cluster/fragment
            mergeBySBR=find(SBR(i,:)> opt.mergeByGeometry.minSBR | SBR(:,i)'> opt.mergeByGeometry.minSBR);
            
            if ~isempty(mergeBySBR)
                %Find the edges that attach to mGrain(i)
                nId_clust=find(i==parentId);
                eAllID_clust=nonzeros(eEdgeId(nId_clust,:));
                
                %We only want to add the minimum number of edges needed for merging, or all the twin
                %type addeds from a relaxed tolerance
                for j=1:length(mergeBySBR)
                    %Find the edges that attach to grains that are set for merging
                    nId = arrlocarr(mergeBySBR(j),parentId);
                    eAllID=nonzeros(eEdgeId(nId,:));

                    %Find the common edges shared between the main cluster mGrain(i) and clusters to
                    %merge
                    toAdd=arrlocarr(eAllID,eAllID_clust);
                    eAllIDAdd=eAllID_clust(toAdd);
                    
                    %If the a twin type boundary(s) exist use those, else use first boundary that
                    %enables the merging.
%                     isTwin=GBTypeRlx(eAllIDAdd)>0;
%                     if any(isTwin)
%                         combineMergeByGeometry(eAllIDAdd(isTwin))=true;
%                     else
%                         combineMergeByGeometry(eAllIDAdd(1))=true;
%                     end
                    if ~isempty(eAllIDAdd)
%                         combineMergeByGeometry(eAllIDAdd(1))=true;
%                         mnId = parentId(eAllIDAdd(1));
                        meId=meEdgeId(parentId(pairs(eAllIDAdd(1),1)),parentId(pairs(eAllIDAdd(1),2)));
                        meGID(meId)=eAllIDAdd(1);
                        mnSBR(mergeBySBR(j))=true;
                        mnType(mergeBySBR(j))=1;
                    end
                end
            end
            
            if (G_geo.Nodes.small(i) || G_geo.Nodes.largeAR(i)) && mnType(i)==0
                %find neighboring merged grains
                indNeighbors1=find(mPairs(:,1)==i);
                indNeighbors2=find( mPairs(:,2)==i);
                nIdmNeighbors=[mPairs(indNeighbors1,2);mPairs(indNeighbors2,1)];
                nIdsmallGrains=find(i==parentId);     

                %find neighbors that contain the orientation
                mergeVote=zeros(length(nIdmNeighbors),1);
                for k=1:length(nIdmNeighbors)
%                     nIdNeighbors=find(nIdmNeighbors(k)==parentId);
                    for j=1:length(nIdsmallGrains)
                        mergeVote(k)=mergeVote(k)+sum(angle(ori(nIdsmallGrains(j)),mFamilyCenters{nIdmNeighbors(k)}) <   opt.mergeByGeometry.FamilyMisTol);
                    end
                end

                %Based on the vote merge neighboring grains
                [valMax,locMax]=max(mergeVote);
                if valMax>0 & mnType(nIdmNeighbors(locMax))==0
                    nIdNeighbors=find(nIdmNeighbors(locMax)==parentId);
                    eAllID_clust=nonzeros(eEdgeId(nIdNeighbors,:));
                    
                    for j=1:length(nIdsmallGrains)
                        %Find the edges that attach to grains that are set for merging
                        eAllID=nonzeros(eEdgeId(nIdsmallGrains(j),:));

                        %Find the common edges shared between the main cluster mGrain(i) and clusters to
                        %merge
                        toAdd=arrlocarr(eAllID,eAllID_clust);
                        eAllIDAdd=eAllID_clust(toAdd);

                        %If the a twin type boundary(s) exist use those, else use first boundary that
                        %enables the merging.
%                         isTwin=GBTypeRlx(eAllIDAdd)>0;
%                         if any(isTwin)
%                             combineMergeByGeometry(find(eAllIDAdd(isTwin))=true;
                        if ~isempty(eAllIDAdd)
%                             combineMergeByGeometry(eAllIDAdd(1))=true;
                            meId=meEdgeId(parentId(pairs(eAllIDAdd(1),1)),parentId(pairs(eAllIDAdd(1),2)));
                            meGID(meId)=eAllIDAdd(1);
                            mnType(i)=2;
                        end
                    end
                end
            end
        end

        %Get the main grains using largeAR and small grains as mask
        G_geo.Edges.combine=meGID~=0;
        G_geo.Edges.eGID=meGID;
        G_geo.Nodes.SBR=mnSBR;
        G_geo.Nodes.eType=mnType;
        G_geo.Nodes.centroids=mGrains.centroid;
        %Remove edges to make reduced graph over the clusters
        G_geo_clust=rmedge(G_geo,G_geo.Edges.EndNodes(~G_geo.Edges.combine,1),...
            G_geo.Edges.EndNodes(~G_geo.Edges.combine,2));


        %Find the large grains
%         G_geo_clust.Nodes.isLarge=~(G_geo_clust.Nodes.largeAR & G_geo_clust.Nodes.small & G_geo_clust.Nodes.SBR);
%         geo_connect = conncomp(G_geo_clust);
%         geo_groups=unique(geo_connect);
%         for i=1:length(geo_groups)
%             %Extract the subgraph
%             mnId=find(geo_connect==i);
%             G_geo_clust_sub = subgraph(G_geo_clust,mnId);
%             
%             %Find any merged grains that include two or more large grains that shouldn't be merged
%             nLarge=sum(G_geo_clust_sub.Nodes.isLarge);
%             if nLarge>1
% %                 %The only problem nodes are those that apply two merging operations on grain
% %                 %We don't allow double merging from shareboundary and small grains, but double
% %                 %merging can occur for 
% %                 nType2=find(G_geo_clust_sub.Nodes.nType==2);
% %                 
% %                 
% %                 meType(meId)=1;
%                 G_geo_clust=rmnode(G_geo_clust,mnId);
%             end
%         end
% 
%         %Plot full graph
%         labelNodes=true;labelEdges=false;plotG=true;legendOn=false;
%         fhandle = plotGraph(grains,mGrains,G_geo_clust,...
%          grains.meanOrientation,1:length(grains),...
%          labelNodes,labelEdges,legendOn,plotG,[]);
        
        %Two large grains should not be combined. To avoid merging to the wrong cluster, the merge
        %process is performed iteratively, and when two large grains are in a cluster, they are both
        %removed. Next iteration the merging direction is determined. 
        combineMergeByGeometry(G_geo.Edges.eGID(G_geo.Edges.combine))=true;
        if any(combineMergeByGeometry)
            
            %Clean the type
            typeTmp=G_frag.Edges.GBTypeRlx(combineMergeByGeometry);
            typeTmp(typeTmp==0)=opt.moriUnknown;
            
            %Store in main arrays
            G_frag.Edges.type(combineMergeByGeometry)=typeTmp;
            G_frag.Edges.combineCleaned(combineMergeByGeometry)=true;

        end
        cnt=cnt+1;
    end %merge inclusion cluster
end

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

        FamilyID{i}=calcClusterAnisotropic(mOri(ngroupId),opt.gfam.alphaFam,opt.gfam.alphaFrag); 
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
end

function [Family] = calcClusterAnisotropic(ori,alphaFam,alphaFrag)

    [cId,fCenters] = calcCluster(ori,'maxAngle',alphaFam,'method','hierarchical');

    %Get the misorientation between family
    while 1
        FamilyMis=round(triu(angle_outer(fCenters,fCenters)./degree),4);
        [r,c]=find(FamilyMis~=0);
 
        if isempty(r) || isempty(c), break; end
        unchangedCnt=0;
        for j=1:size(r,1)
            rid=find(cId==r(j));
            cid=find(cId==c(j));
            if ~isempty(rid) && ~isempty(cid)
                orir=ori(rid);
                oric=ori(cid);

                d = full(abs(dot_outer(orir,oric)));
                dr = full(abs(dot_outer(orir,orir)));
                dc = full(abs(dot_outer(oric,oric)));

                dstart=d;
                cIdstart=cId;
                % progress(0,length(ori));
                while 1

                  % find smallest pair
                  [omega,id] = max(d(:));

                  if omega<cos(alphaFrag), break; end

                  [l,m] = ind2sub(size(d),id);

                  cId(cid(m))=cId(rid(l));

                  d(l,m)=0;
                end
                if all(d(:)==dstart(:)) && all(cId(:)==cIdstart(:))
                    unchangedCnt=unchangedCnt+1;
                end

            else
                unchangedCnt=unchangedCnt+1;
            end
        end
        if unchangedCnt==size(r,1), break; end
    end
    Family=cId;
    
    %relabel the families in case one or more has been fully merged
    cnt=1;
    for i=1:max(Family)
        ind=i==Family;
        if any(ind)
            Family(ind)=cnt;
            cnt=cnt+1;
        end
    end
    
end

function [G_clust,G_frag] = SetTypeForUntwinnedCluster(G_clust,G_frag)
%SetTypeForUntwinnedCluster handles the case were low angle grain
%boundaries are not merged when the fragment graph is made.

    groups = G_clust.Nodes.Group;
    FamilyID = G_clust.Nodes.FamilyID;
    type=G_clust.Nodes.type;
    groupList=unique(groups);
    for i=1:length(groupList)
        group=groupList(i);
        %Find the number of Families
        isGroup=group==groups;
        nFam = max(FamilyID(isGroup));
        if nFam==1
            type(isGroup)=0;
        end
    end
    G_clust.Nodes.type=type;
    G_frag.Nodes.type=type;
end