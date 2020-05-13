function [G_Family,G_clust,G,mGrains] = Cluster(G,grains,opt)
    %Cluster performs the grouping based on edge definitions and user input
    %files. 
    
    % Test if all the files exist and generate if they don't
    flist={'eAddList.txt','eRemoveList.txt','eRemoveCleanup.txt','nRemoveList.txt','nAddList.txt','eRelation.txt','notParent.txt','notTwin.txt'};
    generateInputs(flist); %move to cluster
    
    eAddList=load(flist{1}); %Adds edge and tries to give it a twin relationship, otherwise type twin unknown
    eRemoveList=load(flist{2}); %Remove the edge 
    eRemoveCleanup=load(flist{3}); %Edges determined by cleanFamilyTree
    nRemoveList=load(flist{4}); %Removes all edges associated with node
    nAddList=load(flist{5}); %Tries to add all edges that have twin relationship
        
    %Assign the type using the relaxed mean misorientation tolerance.
    meanTypeRlx=G.Edges.meanTypeRlx;
    G.Edges.type=G.Edges.meanTypeRlx;
    
    %As a starting point, use tight boundary misorientation merging 
    G.Edges.combineCleaned = G.Edges.combine;
    
    %Any edges that aren't a twin should be set to be removed
    %Filter edge by mistol based with the type of meanMistolRelaxed
    G.Edges.combineCleaned(G.Edges.type==0)=false;
    
    %Try adding nodes
    pairs=G.Edges.pairs;
    for i=1:length(nAddList)
        toAdd=any(nAddList(i)==pairs,2);
        G.Edges.combineCleaned(toAdd)=G.Edges.combineRlx(toAdd);
    end
    
    %Try removing nodes
    for i=1:length(nRemoveList)
        [~,epairs]=grains(nRemoveList(i)).neighbors;
        for j=1:size(epairs,1)
            eId=find(all(epairs(j,:)==G.Edges.pairs,2) | all(fliplr(epairs(j,:))==G.Edges.pairs,2));
            isTwin=G.Edges.type(eId)~=0 & G.Edges.type(eId)~=opt.twinUnknown;
            G.Edges.combineCleaned(eId(isTwin))=false;
        end
    end
    
    %Add edges
    eAddId=cell(size(eAddList,1),1);
    for i=1:size(eAddList,1)
        eAddId{i}=find(all(eAddList(i,:)==G.Edges.pairs,2) | all(fliplr(eAddList(i,:))==G.Edges.pairs,2));
    end
    eAddId=cell2mat(eAddId);
    typeAdd=meanTypeRlx(eAddId);
    typeAdd(typeAdd==0)=opt.twinUnknown; %Unknown Type
    G.Edges.type(eAddId)=typeAdd;
    G.Edges.combineCleaned(eAddId)=true;
    
    %Update unknown types in nodes
    typeKnown=unique(G.Edges.pairs([G.Edges.type~=opt.twinUnknown,G.Edges.type~=opt.twinUnknown]));
    G.Nodes.typeUnknown(:)=true;
    G.Nodes.typeUnknown(typeKnown)=false;
  
    %Remove edges
    G.Edges.combineCleaned(eRemoveList)=false; %
    G.Edges.combineCleaned(eRemoveCleanup)=false; %

    
    %Asign each pair to a group of grains
%     for i=1:nMergedGrains
%         nId=find(parentId==i);
%         nGroup(nId)=i;
%         [~,loc]=intersect(pairs(:,1),nId,'rows');
%         eGroup(loc)=i;
%     end
    

    cnt=1;
    ori=grains.meanOrientation;
    [~,ePairs]=grains.neighbors;
    isGNeighbor=zeros(length(grains),length(grains),'logical');
    ind = sub2ind(size(isGNeighbor), ePairs(:,1), ePairs(:,2));
    isGNeighbor(ind)=true;
    
    while cnt<3    
        %Merge grains based on edge combine list
        [mGrains,parentId,mTwinGb,G.Edges.GbInd] = MergeByEdge(G.Edges.pairs,G.Edges.GbInd,G.Edges.combineCleaned,grains);    
        if ~opt.mergeInclusionCluster
            break;
        else
            %Add grains that are mostly internal to a specific cluster to that
            %cluster
            SBR = SharedBoundaryRatio(mGrains);

            %Try adding nodes using any edge that has a twin relationship with
            %Relaxed tolerance else add one edge of unknown type
            pairs=G.Edges.pairs;
            [~,mPairs]=neighbors(mGrains)
            smallmGrains=mGrains.grainSize<200;
%             ePairs=[];
            for i=1:length(mGrains)
                ind=find(SBR(i,:)>0.7);
                if ~isempty(ind)
                    isPair=zeros(size(pairs,1),1,'logical');
%                     nId=grains(i==parentId).id;
                    nId=find(i==parentId);
                    for j=1:length(nId)
                        isPair(nId(j)==pairs(:,1) | nId(j)==pairs(:,2))=true;
                    end
                    pairsInd=find(isPair);
                    pairsSub=pairs(isPair,:);

                    isPairSub=zeros(size(pairsSub,1),1,'logical');
%                     nId=grains(ind==parentId).id;
                    nId=find(ind==parentId);
                    for j=1:length(nId)
                        isPairSub(nId(j)==pairsSub(:,1) | nId(j)==pairsSub(:,2))=true;
                    end

                    toAdd=pairsInd(isPairSub);
                    combine=G.Edges.combineRlx(toAdd);
                    if ~isempty(combine) && ~any(combine)
                       %Then use the first unknown twin relationship to merge into
                       %cluster. Note this addition will go away if nAddList doesn't
                       %contain the same node in futures calls to cluster or the 
                       %first edge is in edge remove list.
                       combine(1)=true; 
                       G.Edges.type(toAdd(1))=opt.twinUnknown;
                    end
                    G.Edges.combineCleaned(toAdd)=combine;
                end
                if smallmGrains(i)
                    %find neighboring merged grains
                    indNeighbors=find(mPairs(:,1)==i);
                    nIdmNeighbors=mPairs(indNeighbors,2);
                    nIdsmallGrains=find(i==parentId);     
                    
                    %find neighbors that contain the orientation
                    mergeVote=zeros(length(nIdmNeighbors),1);
                    for k=1:length(nIdmNeighbors)
                        nIdNeighbors=find(nIdmNeighbors(k)==parentId);
                        for j=1:length(nIdsmallGrains)
                            mergeVote(k)=mergeVote(k)+sum(angle(ori(nIdsmallGrains(j)),ori(nIdNeighbors)) < opt.grain_recon.seg_angle_grouped);
                        end
                    end
                    
                    %Based on the vote merge neighboring grains
                    [valMax,locMax]=max(mergeVote);
%                     fprintf('grain number %d of %d\n',i,length(mGrains))
                    if valMax>0
                        nIdNeighbors=find(nIdmNeighbors(locMax)==parentId);
%                         gBId2=grains(nIdNeighbors).boundary.grainId;
                        for k=1:length(nIdsmallGrains)
%                             [row,~]=find(nIdsmallGrains(k)==gBId2);
%                             uniqueNeighbors=unique(gBId2(row,:));
%                             uniqueNeighbors(uniqueNeighbors==nIdsmallGrains(k))=[];
                            uniqueNeighbors=nIdNeighbors(isGNeighbor(nIdsmallGrains(k),nIdNeighbors));
                            eAddList=vertcat(eAddList,[uniqueNeighbors,nIdsmallGrains(k)*ones(length(uniqueNeighbors),1)]);
                        end
                    end
                end
                
            end
            
            if ~isempty(eAddList)
                %Add edges
                eAddId=cell(size(eAddList,1),1);
                for i=1:size(eAddList,1)
                    eAddId{i}=find(all(eAddList(i,:)==G.Edges.pairs,2) | all(fliplr(eAddList(i,:))==G.Edges.pairs,2));
                end
                eAddId=cell2mat(eAddId);
                typeAdd=meanTypeRlx(eAddId);
                typeAdd(typeAdd==0)=opt.twinUnknown; %Unknown Type
                G.Edges.type(eAddId)=typeAdd;
                
                G.Edges.combineCleaned(eAddId)=true;

                %Update unknown types in nodes
                typeKnown=unique(G.Edges.pairs([G.Edges.type~=opt.twinUnknown,G.Edges.type~=opt.twinUnknown]));
                G.Nodes.typeUnknown(:)=true;
                G.Nodes.typeUnknown(typeKnown)=false;

                %Remove edges
                G.Edges.combineCleaned(eRemoveList)=false; 
                G.Edges.combineCleaned(eRemoveCleanup)=false; 
            end
        end
        cnt=cnt+1;
    end %merge inclusion cluster
    
   [mGrains,parentId,mTwinGb,G.Edges.GbInd] = MergeByEdge(G.Edges.pairs,G.Edges.GbInd,G.Edges.combineCleaned,grains);    

    
    %Try removing nodes so no internal grains are added that we specify to
    %remove
    for i=1:length(nRemoveList)
        [~,epairs]=grains(nRemoveList(i)).neighbors;
        for j=1:size(epairs,1)
            eId=find(all(epairs(j,:)==G.Edges.pairs,2) | all(fliplr(epairs(j,:))==G.Edges.pairs,2));
            isTwin=G.Edges.type(eId)~=0 & G.Edges.type(eId)~=opt.twinUnknown;
            G.Edges.combineCleaned(eId(isTwin))=false;
        end
    end

    %Remove edges similarly to ensure nothing is added that we said not to
    %add
    G.Edges.combineCleaned(eRemoveList)=false; %
    G.Edges.combineCleaned(eRemoveCleanup)=false;
        
    %Remove edges to make reduced graph over the clusters
    G_clust=rmedge(G,G.Edges.pairs(~G.Edges.combineCleaned,1),...
        G.Edges.pairs(~G.Edges.combineCleaned,2));
    
    %Assign group Ids using mGrains id
    [G_clust.Nodes.Group,G_clust.Edges.Group] = GroupsFromGraph(G_clust.Edges.pairs,length(grains),length(mGrains),parentId);
    
    %Transfer still valid Family and family vote info to G_clust.. build
    %list of groups that need to be computed
    [groupList] = GetGroupsToCompute(G_clust,G);
    
    %Update G with groups in G_clust
    G.Nodes.Group(G_clust.Nodes.Id)=G_clust.Nodes.Group;
    G.Edges.Group(G.Edges.combineCleaned)=G_clust.Edges.Group;
    
    %Ensure no zero types (DS: This shouldn't happen if logic is airtight.. which means something is off somewhere.)
    G_clust.Edges.type(G_clust.Edges.type==0)=opt.twinUnknown;
    
    %Update unknown nodes
    nId_edges=unique(G_clust.Edges.pairs([G_clust.Edges.type~=opt.twinUnknown,G_clust.Edges.type~=opt.twinUnknown]));
    G_clust.Nodes.typeUnknown(:)=true; 
    G_clust.Nodes.typeUnknown(nId_edges)=false;

    % Identify Families in grain clusters 
    [G_clust,G]= AssignFamilyIDs(G_clust,G,grains,mGrains,opt);


    %Create the minimum set of pairs that should be tested for twin
    %relationships for each family group size
    pairIndexes=cell(max(G_clust.Nodes.FamilyID),1);
    pairSize=zeros(max(G_clust.Nodes.FamilyID),1);
    for numFamily=1:max(G_clust.Nodes.FamilyID)
        [row,col]=find(triu(ones(numFamily))-eye(numFamily));
        pairIndexes{numFamily}=[row,col];
        pairSize(numFamily)=length(row);
    end
    
    % Construct Family twin relationship graph 
    numFamilypGroup=zeros(length(mGrains),1)
    for i=1:length(mGrains)
        group=mGrains.id(i);
        ngroupId = find((group==G_clust.Nodes.Group)==true);
        numFamilypGroup(i) = double(max(G_clust.Nodes.FamilyID(ngroupId)));
    end
    
    %Initialize nodes
    numEdges=sum(pairSize(numFamilypGroup));
    s=zeros(numEdges,1);
    t=zeros(numEdges,1);
    Group=zeros(numEdges,1);
    familyPair=zeros(numEdges,2);

    %Initialize graph
    ind=1;
    lastPairId=0;
    for i=1:length(mGrains)
        if pairSize(numFamilypGroup(i))~=0
            pInd=pairIndexes{numFamilypGroup(i)};
            indx=ind:ind-1+size(pInd,1);
            s(indx)=pInd(:,1)+lastPairId;
            t(indx)=pInd(:,2)+lastPairId; %add lastPairId so no conflicts
            Group(indx)=i;
            familyPair(indx,:)=pInd;
            ind=indx(end)+1;
            lastPairId=lastPairId+numFamilypGroup(i);
        end
    end
    G_Family=graph(s,t);
    G_Family.Edges.pairs=[s,t];
    G_Family.Edges.Group=Group;
    G_Family.Edges.familyPair=familyPair;
    
    % Construct node Family and Group
    [~,ind]=unique(G_Family.Edges.pairs)
    G_Family.Nodes.Family=familyPair(ind);
    nGroup(s)=Group;
    nGroup(t)=Group;
    G_Family.Nodes.Group=nGroup';
    
    % Construct Family twin relationship graph 
    oriAll=grains.meanOrientation ;
    areaAll=grains.area;
    maxNumFamilies=max(numFamilypGroup)
    ori1 = orientation.byEuler(zeros(numEdges,1),...
        zeros(numEdges,1),zeros(numEdges,1),...
            'ZYZ',opt.CS{2})
    ori2=ori1;
    eFType=zeros(numEdges,1);
    for i=1:length(mGrains)
        ngroupId = find(i==G_clust.Nodes.Group);
        egroupId = find(i==G.Edges.Group); %converts logical arrays to indices
        egroupFId = find(i==G_Family.Edges.Group);
        nFamily = G.Nodes.FamilyID(ngroupId);     
        eFamily = G.Edges.FamilyID(egroupId,:);
        eType = G.Edges.type(egroupId);
        familyPair=G_Family.Edges.familyPair(egroupFId,:);

        %transfer edge based types to family types
        for j=1:length(egroupFId)
            ind=all(eFamily==familyPair(j,:) | fliplr(eFamily)==familyPair(j,:),2);
            eTypeLoc=eType(ind);
            eTypeLoc=eTypeLoc(eTypeLoc<opt.twinUnknown & eTypeLoc>0);
            if ~isempty(eTypeLoc)
                eFType(egroupFId(j))=mode(eTypeLoc);
            end
        end
        
        for j=1:max(nFamily)
            oriloop=oriAll(ngroupId(nFamily==j));
            arealoop=areaAll(ngroupId(nFamily==j));
            ori1(egroupFId(familyPair(:,1)==j)) = mean(oriloop,'weights',arealoop);            
            ori2(egroupFId(familyPair(:,2)==j)) = mean(oriloop,'weights',arealoop);
        end 
    end
    
    mori=inv(ori1).*ori2; 
    oriFamily=[ori1,ori2]
    tol=zeros(opt.nTwin,1);
    tolRlx=zeros(opt.nTwin,1);
    for i=1:opt.nTwin
        tol(i)=opt.twin{i}.tol.misMean;
        tolRlx(i)=opt.twin{i}.tol.misMeanRlx;
    end
    meanType=zeros(numEdges,1);
    meanTypeRlx=zeros(numEdges,1);
    [~,G_Family.Edges.meanType] = TestTwinRelationship(mori,tol,opt,eFType);
    [~,G_Family.Edges.meanTypeRlx] = TestTwinRelationship(mori,tolRlx,opt,eFType);
    
    %Compute Schmid info for twin/parents in clustered grains
    %This computes Schmid factor for twin/parent identification and 
    %is stored at both edge and node level
    [G_Family,G_clust]=GetSchmidRelative(G_Family,G_clust,oriFamily,G_Family.Edges.meanTypeRlx,grains,mGrains,opt)
%     [G_clust,G]= GetSchmidRelative(G_clust,G,grains,mGrains,opt);
    
    %Perform family votes
    [G_Family,G_clust] = FamilyVotes(G_Family,G_clust,unique(G_clust.Nodes.Group),grains,mGrains,opt);    

    %Plot the results
    if opt.plot.do
        %Remove grains that haven't twinned (better plot performance)
        if opt.plot.ClusterOnly
            tokeep=zeros(length(G_clust.Nodes.Id),1,'logical');
            tokeep(unique(G_clust.Edges.pairs))=true;
            G_clust_clean=rmnode(G_clust,find(~tokeep));
        else
            G_clust_clean=G_clust;
        end
        
        %Plot edge labeled graph
        labelNodes=false;labelEdges=opt.plot.labelEdges;plotG=true;legendOn=false;
        fhandle = plotGraph(grains,mGrains,G_clust_clean,...
            grains.meanOrientation,G_clust_clean.Nodes.Id,...
            labelNodes,labelEdges,legendOn,plotG,[]);
        
        %Plot node labeled graph
        labelNodes=true;labelEdges=false;plotG=false;legendOn=false;
        fhandle = plotGraph(grains,mGrains,G_clust_clean,...
            grains.meanOrientation,G_clust_clean.Nodes.Id,...
            labelNodes,labelEdges,legendOn,plotG,[]);
        
    end
end
function [groupList] = GetGroupsToCompute(G_clust,G)
    %Get group ids - optimized 12/24/19
    groupCompute=zeros(max(G_clust.Nodes.Group),1,'logical');
    for i=1:max(G_clust.Nodes.Group)
        %Extract group to local variables for readability
        eID=G_clust.Edges.GlobalID(i==G_clust.Edges.Group);
        nID = find((i==G_clust.Nodes.Group));
        G_nID = find(G.Nodes.Group(nID(1))==G.Nodes.Group);
        G_eId = find(G.Nodes.Group(nID(1))==G.Edges.Group);
               
        if length(G_nID) ~= length(nID) || length(G_eId) ~= length(eID) 
            groupCompute(i)=true;
        elseif ~all(nID==G_nID) || ~all(eID==G_eId) %Because sorted
            groupCompute(i)=true;           
        end
    end
    groupList=find(groupCompute);
end
