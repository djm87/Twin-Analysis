function [G_clust,G,mGrains] = ClusterGraph(G,grains,opt)
    %Cluster performs the grouping based on edge definitions and user input
    %files. 
    
    % Test if all the neccessary files exist and generate if they don't
    flist={'eAddList.txt','eRemoveList.txt','nRemoveList.txt','nAddList.txt','eRelation.txt','notParent.txt','notTwin.txt'};
    generateInputs(flist); %move to cluster
    
    % Load content of files
    eAddList=load(flist{1}); %Adds edge and tries to give it a twin relationship, otherwise type twin unknown
    eAddList(eAddList>G.numedges)=[];
    eRemoveList=load(flist{2}); %Remove the edge 
    eRemoveList(eRemoveList>G.numedges)=[];
    nRemoveList=load(flist{3}); %Removes all edges associated with node
    nAddList=load(flist{4}); %Tries to add all edges that have twin relationship
        
    %As a starting point, use tight boundary misorientation. Both
    %combine and GBType are not changed after initialization
    G.Edges.combineCleaned = G.Edges.combine;
    G.Edges.type=G.Edges.typeIni;
    
    %To simplify edge id extraction from a node pair
    pairs=G.Edges.pairs; 
    maxEdgeId=size(pairs,1);
    AllPairs=vertcat(pairs,fliplr(pairs));
    AllPairId=double(vertcat(G.Edges.GlobalID,G.Edges.GlobalID));
    eEdgeId=sparse(AllPairs(:,1),AllPairs(:,2),AllPairId,maxEdgeId,maxEdgeId);
    
    %Try adding all edges attached to a user specified node using a relaxed tolerance
    eToAdd=nonzeros(eEdgeId(nAddList,:));
    eToAdd=eToAdd(G.Edges.combineRlx(eToAdd));
    G.Edges.combineCleaned(eToAdd)=true;
    G.Edges.type(eToAdd)=G.Edges.GBTypeRlx(eToAdd);
    
    %Add user specified edge
    eToAdd=zeros(size(eAddList,1),1);
    for i=1:size(eAddList,1)
        tmp=nonzeros(eEdgeId(eAddList(i,1),eAddList(i,2)));
        if ~isempty(tmp)
            eToAdd(i)=tmp;
        end
    end
    eToAdd(eToAdd==0)=[];
    eTypeAdd=G.Edges.meanType(eToAdd);
    eTypeAdd(eTypeAdd==0)=opt.twinUnknown; %unknown type
    G.Edges.type(eToAdd)=eTypeAdd;
    G.Edges.combineCleaned(eToAdd)=true;
    
    %Update unknown types in nodes
    nId_edges=unique(G.Edges.pairs([G.Edges.type~=opt.twinUnknown,G.Edges.type~=opt.twinUnknown]));
    G.Nodes.typeUnknown(:)=true;
    G.Nodes.typeUnknown(nId_edges)=false;
    
    %Removing all edges attached to a user specified node
    toRemove=nonzeros(eEdgeId(nRemoveList,:));
    G.Edges.combineCleaned(toRemove)=false;
    
    %remove user specified edge
    eRemoveList(eRemoveList>G.numedges)=[];
    G.Edges.combineCleaned(eRemoveList)=false;
    
    %Add back the geometry based merging when there is some geometry to
    %test it on
    if opt.mergeByGeometry.mergeCluster
        G = MergeByGeometry(G,eEdgeId,grains,opt)
    end
    
    %Update unknown types in nodes
    nId_edges=unique(G.Edges.pairs([G.Edges.type~=opt.twinUnknown,G.Edges.type~=opt.twinUnknown]));
    G.Nodes.typeUnknown(:)=true;
    G.Nodes.typeUnknown(nId_edges)=false;
    
    %Removing all edges attached to a user specified node
    toRemove=nonzeros(eEdgeId(nRemoveList,:));
    G.Edges.combineCleaned(toRemove)=false;
    
    %remove user specified edge
    G.Edges.combineCleaned(eRemoveList)=false;
    
    %Remove edges to make reduced graph over the clusters
    G_clust=rmedge(G,G.Edges.pairs(~G.Edges.combineCleaned,1),...
        G.Edges.pairs(~G.Edges.combineCleaned,2));
    
    %Ensure no zero types (This should always be true unless there is an error in logic).
    assert(all(G_clust.Edges.type~=0),'There is an error in type definition in clusterGraph')
    
    %Build the merged grains from the edge graph
   [mGrains,parentId,mTwinGb,G.Edges.GbInd] = MergeByEdge(G.Edges.pairs,G.Edges.GbInd,G.Edges.combineCleaned,grains);    
    
    %Assign group Ids using mGrains id
    [G_clust.Nodes.Group,G_clust.Edges.Group] = GroupsFromGraph(G_clust.Edges.pairs,length(grains),length(mGrains),parentId);
    
    %Transfer still valid Family and family vote info to G_clust.. build
    %list of groups that need to be computed
    [G_clust,G] = GetGroupsToCompute(G_clust,G);
    
    %Update G with groups in G_clust
    G.Nodes.Group(G_clust.Nodes.Id)=G_clust.Nodes.Group;
    G.Edges.Group(:)=0;
    G.Edges.Group(G.Edges.combineCleaned)=G_clust.Edges.Group;
    
    %Update unknown nodes
    nId_edges=unique(G_clust.Edges.pairs([G_clust.Edges.type~=opt.twinUnknown,G_clust.Edges.type~=opt.twinUnknown]));
    G_clust.Nodes.typeUnknown(:)=true; 
    G_clust.Nodes.typeUnknown(nId_edges)=false;

    % Identify Families in grain clusters 
    [G_clust,G]= AssignFamilyIDs(G_clust,G,grains,mGrains,opt);
    
    % Enforce clusters with single Family to have type = 0
    % This is enforced because types are saved onced the Family has been computed and comes back
    % clean
    [G_clust,G] = SetTypeForUntwinnedCluster(G_clust,G);
end
function [G_clust,G] = GetGroupsToCompute(G_clust,G)
    %Get group ids to process
    %Notes: Initially G.Nodes.Group is true and all clusters will be
    %triggered for computation. Clusters haven't been computed until
    %G_clust comes out of the Family graph routine. Similarly computeFamily
    %isn't complete until a cleanFamilyTree returns exflag=0 for a cluster.
    %The point of this routine is to limit recomputation, but also allow
    %for changes to Family graph to happen once based off some user
    %interaction and stay as long as the clusters don't have grains added
    %or removed. G_clust is anchored in G as long as the EBSD segmentation
    %doesn't change, G_Family is not anchored in G_clust in general since
    %the grain grouping ends up being an important aspect of creating the
    %family tree. IsNewGroup and computeFamily are seperate since one ends
    %with creating the Family graph and the other ends when a valid family
    %tree is produced.
%     groupCompute=zeros(max(G_clust.Nodes.Group),1,'logical');

    for i=1:max(G_clust.Nodes.Group)
        %For a given G_clust cluster get one grain
        eID=G_clust.Edges.GlobalID(i==G_clust.Edges.Group);
        nID = find((i==G_clust.Nodes.Group));
        
        %Find the stored group in G that contains that grain
        G_nID = find(G.Nodes.Group(nID(1))==G.Nodes.Group);
        G_eId = find(G.Nodes.Group(nID(1))==G.Edges.Group);
        
        %If the store node ids match the node ids of G_clust, then the
        %cluster has not changed
        if length(G_nID) ~= length(nID) || length(G_eId) ~= length(eID) 
%             groupCompute(i)=true;
            G_clust.Nodes.isNewGroup(nID)=true;
            G_clust.Nodes.computeFamily(nID)=true;
        elseif ~all(nID==G_nID) || ~all(eID==G_eId) %Because sorted
%             groupCompute(i)=true; 
            G_clust.Nodes.isNewGroup(nID)=true;
            G_clust.Nodes.computeFamily(nID)=true;
        end
    end
%     groupList=find(groupCompute);
    %Store in persistent graph
    G.Nodes.isNewGroup(G_clust.Nodes.Id)=G_clust.Nodes.isNewGroup;
    G.Nodes.computeFamily(G_clust.Nodes.Id)=G_clust.Nodes.computeFamily;

    %Reset the Family information
    G_clust.Nodes.Generation(G.Nodes.isNewGroup)=-1;
    G.Nodes.Generation=G_clust.Nodes.Generation;
    G_clust.Nodes.type(G.Nodes.isNewGroup)=0;
    G.Nodes.type=G_clust.Nodes.type;
    G_clust.Nodes.ParentFamilyId(G.Nodes.isNewGroup)=0;
    G.Nodes.ParentFamilyId =G_clust.Nodes.ParentFamilyId;
    G_clust.Nodes.SchmidVariant(G.Nodes.isNewGroup)=0;
    G.Nodes.SchmidVariant=G_clust.Nodes.SchmidVariant;
    G_clust.Nodes.K1NAng(G.Nodes.isNewGroup)=0;
    G.Nodes.K1NAng=G_clust.Nodes.K1NAng;
    
end
function [G] = MergeByGeometry(G,eEdgeId,grains,opt)
       
    %Initialize variables
    ori=grains.meanOrientation;
%     isGNeighbor=zeros(length(grains),length(grains),'logical');
%     ind = sub2ind(size(isGNeighbor), ePairs(:,1), ePairs(:,2));
%     isGNeighbor(ind)=true;

%     pairs=G.Edges.pairs; 
%     maxEdgeId=size(pairs,1);
%     AllPairs=vertcat(pairs,fliplr(pairs));
%     AllPairId=double(vertcat(G.Edges.GlobalID,G.Edges.GlobalID));
%     eEdgeId=sparse(AllPairs(:,1),AllPairs(:,2),AllPairId,maxEdgeId,maxEdgeId);
    
    ePairs=G.Edges.pairs;
    cnt=1;
    while cnt<=opt.mergeByGeometry.itter    
        %Merge grains based on edge combineCleaned list
        [mGrains,parentId,mTwinGb,G.Edges.GbInd] = MergeByEdge(G.Edges.pairs,G.Edges.GbInd,G.Edges.combineCleaned,grains);    

        %Plot full graph
%         labelNodes=true;labelEdges=false;plotG=false;legendOn=false;
%         fhandle = plotGraph(grains,mGrains,G,...
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
        GBTypeRlx=G.Edges.GBTypeRlx;
        pairs=G.Edges.pairs;
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
                    nIdNeighbors=find(nIdmNeighbors(k)==parentId);
                    for j=1:length(nIdsmallGrains)
                        mergeVote(k)=mergeVote(k)+sum(angle(ori(nIdsmallGrains(j)),ori(nIdNeighbors)) <  opt.mergeByGeometry.FamilyMisTol);
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
            typeTmp=G.Edges.GBTypeRlx(combineMergeByGeometry);
            typeTmp(typeTmp==0)=opt.twinUnknown;
            
            %Store in main arrays
            G.Edges.type(combineMergeByGeometry)=typeTmp;
            G.Edges.combineCleaned(combineMergeByGeometry)=true;

        end
        cnt=cnt+1;
    end %merge inclusion cluster

   
end

function [G_clust,G] = SetTypeForUntwinnedCluster(G_clust,G)
    
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
    G.Nodes.type=type;
end