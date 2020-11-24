function [G_family,G_clust,G_frag] = FamilyGraph(G_family,G_clust,G_frag,grains,opt,recompute,computemGrainId)
    %FamilyGraph creates a graph for each merged grain in the cluster graph
    %and considers all combinations of relationships between families. The
    %nodes are a families and the edges a relationship between families.
    %Each relationship is tested for being a twin, the parent is determined
    %based on the voting scheme of Marshel et al. This routine essentially
    %computes the information needed to solve the family tree.

    %Get groups to compute from G_family,recompute,and computemGrainId
    [groupList,groups,nGroups,toCompute] =...
        getGroupList(G_family,G_clust,recompute,computemGrainId);
    if isempty(groupList)
        fprintf('There is no cluster to compute family graph for!\n')
        return
    end

    %Get the family ID for each node of in clusters to compute
    FamilyID=G_clust.Nodes.FamilyID;
    maxFamilyID=max(FamilyID(toCompute));
    
    %Create the set of pairs that should be tested for twins
    %relationships for each family group size
    [pairIndexes,pairSize] = getFamilyCombinations(maxFamilyID);
    
    % Get the number of families for each group 
    nFamilyGroup = getNFamilyGroup(groups,groupList,FamilyID,nGroups);
    
    % Initialize graph
    [G_family,nEdges] = InitializeFamilyGraph(G_family,pairSize,pairIndexes,nFamilyGroup,groupList,nGroups);
    if nEdges==0
       fprintf('No cluster needs family graph computed!\n')
       return 
    end
    
    % Construct Family twin relationship graph 
    [G_family,oriFamily,mori,eFType] = getOriFamily(G_family,G_clust,grains,groupList,nGroups,nEdges,opt);

    % Test Family pairs for twin relationship
    [~,G_family.Edges.meanType] = TestTwinRelationship(mori,opt.gfam.phi,opt,eFType);
    
    %Compute Schmid info for twin/parents in families
    [G_family,edgeList] = GetSchmidRelative(G_family,groupList,oriFamily,G_family.Edges.meanType,opt);
    
    %Perform family votes
    [G_family] = FamilyVotes(G_family,G_clust,groupList,grains,opt);    

    %Get the parent relationships based on votes
    [G_family,G_clust] = getParent(groupList,G_family,G_clust) 
    
    %Store the schmid factor in graphs G_frag and G_clust for plotting
    [G_clust,G_frag] = StoreSchmid(G_family,G_clust,G_frag,groupList,opt);
    
    %Update isNewGroup so computation doesn't happen again.
    G_clust.Nodes.isNewGroup(:)=false;
    G_frag.Nodes.isNewGroup=G_clust.Nodes.isNewGroup;

end
function [groupList,groups,nGroups,toCompute] = getGroupList(G_family,G_clust,recompute,computemGrainId)
    %Returns: 
    %groupList (nCluster x 1) list of group numbers to compute 
    %groups (nNode x 1) list where each node has a group number
    %nGroups length(groupList) is the number of groups that need computing
    %toCompute (nNode x 1) logical array for nodes with groups that will be
    %computed
   
    nClusters=length(G_family);
    groups=G_clust.Nodes.Group;

    if recompute 
        %Recompute all clusters
        toCompute=ones(length(groups),1,'logical')
    else
        toCompute=zeros(length(groups),1,'logical');
        for i=1:nClusters
            if isempty(G_family{i})
                toCompute(i==groups)=true;
            end
        end
        ind=intersect_wrepeat(computemGrainId,G_clust.Nodes.Group);
        toCompute(ind)=true;
    end
    
    groupList=unique(G_clust.Nodes.Group(toCompute));
    nGroups=length(groupList);
end
function [G_family,G_clust] = getParent(groupList,G_family,G_clust) 
    G_family.Edges.Parent = zeros(size(G_family.Edges.pairs),'logical');    
    
    for i=1:length(groupList)
        group=groupList(i);
        %load edge and node properties for Family graph
        G_Family_sub = subgraph(G_family,find(group==G_family.Nodes.Group));
        nFamily = G_Family_sub.Nodes.Family;
        eType = G_Family_sub.Edges.meanType;
        if ~isempty(eType)
            eVote = G_Family_sub.Edges.Vote;
            eFamily = G_Family_sub.Edges.FamilyID;
            eGlobalID = G_Family_sub.Edges.eGlobalID;

            %array of size(ePairs) of the family correlated with the edge Type 
            FamilyRelationList = FamilyRelationships(nFamily,eType,eFamily); 

            %eVotesSummed are the votes for each for a family and a type
            [eVotesSummed,parentID] = buildeVotesSummed(nFamily,eType,eVote,FamilyRelationList);

            %Determine the parent by vote
            Parent = parentByVote(eVotesSummed,eType,eFamily);

            %Store the parent
            G_family.Edges.Parent(eGlobalID,:)=Parent; %Need to add back in
        else
            G_clust.Nodes.computeFamily(group==G_clust.Nodes.Group)=false;
        end
    end
    
    %update the directions of the graph and remove non-twin edges
    try
        toRemove=find(all(G_family.Edges.Parent==0,2));
        G_family=rmedge(G_family,toRemove);
        G_family.Edges.eGlobalID=[1:numedges(G_family)]';
        toFlip=G_family.Edges.Parent(:,2);
        G_family.Edges.pairs(toFlip,:)=fliplr(G_family.Edges.pairs(toFlip,:));
        G_family.Edges.FamilyID(toFlip,:)=fliplr(G_family.Edges.FamilyID(toFlip,:));
        G_family.Edges.eRelationship(toFlip,:)=fliplr(G_family.Edges.eRelationship(toFlip,:));
        G_family.Edges.sigma13(toFlip,:)=fliplr(G_family.Edges.sigma13(toFlip,:));
        G_family.Edges.FRArea(toFlip,:)=fliplr(G_family.Edges.FRArea(toFlip,:));
        G_family.Edges.FRgB(toFlip,:)=fliplr(G_family.Edges.FRgB(toFlip,:));
        G_family.Edges.FREffSF(toFlip,:)=fliplr(G_family.Edges.FREffSF(toFlip,:));
        G_family.Edges.Vote(toFlip,:)=fliplr(G_family.Edges.Vote(toFlip,:));
        G_family.Edges.Parent(toFlip,:)=fliplr(G_family.Edges.Parent(toFlip,:));
        G_family=flipedge(G_family,find(toFlip));
    catch
    end
end
function FamilyRelationList = FamilyRelationships(nFamily,eType,eFamily) 
    %Returns max(nFamily) x max(eType) cell array containing a logical
    %array of size(ePairs) of the family correlated with the edge Type 
    FamilyRelationList = cell(max(nFamily),max(eType));
    for j = 1:max(nFamily)
        for k = 1:max(eType)
            FamilyInPair = ismember(eFamily(:,:),j);
            [r,c]=find(FamilyInPair);
            FamilyInPair(FamilyInPair) = eType(r)==k;
            FamilyRelationList{j,k} = FamilyInPair;
        end
    end
end
function [eVotesSummed,parentID] = buildeVotesSummed(nFamily,eType,eVote,FamilyRelationList)
    eVotesSummed = zeros(max(nFamily),max(eType));
    for j = 1:max(nFamily)
        for k = 1:max(eType)
            eVotesSummed(j,k) = sum(eVote(FamilyRelationList{j,k}(:,:)));
        end
    end
    [~,parentID] = sort(eVotesSummed,'descend'); %sorts each column
end
function Parent = parentByVote(eVotesSummed,eType,eFamily)
    %Initialize parent
    Parent = zeros(size(eFamily,1),2,'logical');
    for j = 1:size(eFamily,1)
        n1=eFamily(j,1);
        n2=eFamily(j,2);
        if eType(j)~=0
            v1=eVotesSummed(n1,eType(j));
            v2=eVotesSummed(n2,eType(j));
            if v1>v2
               Parent(j,1)=true;
            else
               Parent(j,2)=true; 
            end
        end
    end
end
function [pairIndexes,pairSize] = getFamilyCombinations(maxFamilyID)
    pairIndexes=cell(maxFamilyID,1);
    pairSize=zeros(maxFamilyID,1);
    for nFamily=1:maxFamilyID
        [row,col]=find(triu(ones(nFamily))-eye(nFamily));
        pairIndexes{nFamily}=sortrows([row,col]);
        pairSize(nFamily)=length(row);
    end
end
function [nFamilyGroup] = getNFamilyGroup(groups,groupList,FamilyID,nGroups)
    nFamilyGroup=zeros(nGroups,1);
    for i=1:nGroups
        nFamilyGroup(i) = double(max(FamilyID(groupList(i)==groups)));
    end
end
 function [G_family,nEdges] = InitializeFamilyGraph(G_family,pairSize,pairIndexes,nFamilyGroup,groupList,nGroups)

    %rewrite this section 
 
    %Initialize arrays
    nEdges=sum(pairSize(nFamilyGroup));
    s=zeros(nEdges,1);
    t=zeros(nEdges,1);
    Group=zeros(nEdges,1);
    FamilyID=zeros(nEdges,2);

    %build graph based on possible family pairs
    ind=1;
    lastPairId=0;
    for i=1:nGroups
        if pairSize(nFamilyGroup(i))~=0
            pInd=pairIndexes{nFamilyGroup(i)};
            indx=ind:ind-1+size(pInd,1);
            s(indx)=pInd(:,1)+lastPairId;
            t(indx)=pInd(:,2)+lastPairId; %add lastPairId so no conflicts
            Group(indx)=groupList(i);
            FamilyID(indx,:)=pInd;
            ind=indx(end)+1;
            lastPairId=lastPairId+nFamilyGroup(i);
        end
    end
    if isempty(s)
        G_family=[];
    else
        G_family=digraph(s,t);
        G_family.Edges.pairs=[s,t];
        G_family.Edges.Group=Group;
        G_family.Edges.FamilyID=FamilyID;
        G_family.Edges.eRemove=zeros(nEdges,1,'logical');
        G_family.Edges.eIsParentAll=zeros(nEdges,1,'logical');
        G_family.Edges.eNotParentAll=zeros(nEdges,1,'logical');
        G_family.Edges.eNotParent=zeros(nEdges,1,'logical');
        G_family.Edges.eRelationship=zeros(nEdges,2,'int8');
        G_family.Edges.eGlobalID=[1:nEdges]';
        G_family.Nodes.GlobalID=[1:numnodes(G_family)]';
        G_family.Nodes.root=zeros(numnodes(G_family),1,'logical');

        % Construct node Family and Group
        [~,ind]=unique(G_family.Edges.pairs);
        G_family.Nodes.Family=FamilyID(ind);
        nGroup(s)=Group;
        nGroup(t)=Group;
        G_family.Nodes.Group=nGroup';
    end
end
function [G_family,oriFamily,mori,eFType] = getOriFamily(G_family,G_clust,grains,groupList,nGroups,nEdges,opt)
    oriAll=grains.meanOrientation ;
    areaAll=grains.area;
    ori1 = orientation.byEuler(zeros(nEdges,1),...
        zeros(nEdges,1),zeros(nEdges,1),...
            'ZYZ',opt.CS{2});
    ori2=ori1;
    eFType=zeros(nEdges,1);
    for i=1:nGroups
        group=groupList(i);
        ngroupId = find(group==G_clust.Nodes.Group);
        egroupId = find(group==G_clust.Edges.Group); %converts logical arrays to indices
        egroupFId = find(group==G_family.Edges.Group);
        nFamily = G_clust.Nodes.FamilyID(ngroupId);     
        eFamily = G_clust.Edges.FamilyID(egroupId,:);
        eType = G_clust.Edges.type(egroupId);
        FamilyID=G_family.Edges.FamilyID(egroupFId,:);

        %transfer edge based types to family types
        for j=1:length(egroupFId)
            ind=all(eFamily==FamilyID(j,:) | fliplr(eFamily)==FamilyID(j,:),2);
            eTypeLoc=eType(ind);
            eTypeLoc=eTypeLoc(eTypeLoc<opt.twinUnknown & eTypeLoc>0);
            if ~isempty(eTypeLoc)
                eFType(egroupFId(j))=mode(eTypeLoc);
            end
        end
        
        for j=1:max(nFamily)
            oriloop=oriAll(ngroupId(nFamily==j));
            arealoop=areaAll(ngroupId(nFamily==j));
            ori1(egroupFId(FamilyID(:,1)==j)) = mean(oriloop,'weights',arealoop);            
            ori2(egroupFId(FamilyID(:,2)==j)) = mean(oriloop,'weights',arealoop);
        end 
    end
    %Put the mean orientations in the graph for schmid computation
    [~,ind]=unique(reshape(G_family.Edges.pairs,[numel(G_family.Edges.pairs),1]));
    oriList=[ori1;ori2];
    G_family.Nodes.meanOri=oriList(ind);
    
    mori=inv(ori1).*ori2; 
    oriFamily=[ori1,ori2];
end
function [G_clust,G_frag] = StoreSchmid(G_family,G_clust,G_frag,groupList,opt)
    
    %First store new info in G_clust - EffSF
    nMori=opt.nMori;
    for i=1:length(groupList) 
        group=groupList(i);
        eFID = find(group==G_family.Edges.Group);
        if ~isempty(eFID)
%             nFID = find(group==G_family.Nodes.Group);
%             eCID = find(group==G_clust.Edges.Group);
            nCID = find(group==G_clust.Nodes.Group);            
            nID = find(group==G_frag.Nodes.Group);

            eFEffSF = G_family.Edges.EffSF(eFID,:);
            eFType = G_family.Edges.meanType(eFID);
            eFFamilyID=G_family.Edges.FamilyID(eFID,:);
            nCFamilyID = G_clust.Nodes.FamilyID(nCID);
            nFamilyID = G_frag.Nodes.FamilyID(nID);
            
            nCEffSF=zeros(length(nCID),nMori);
            nEffSF=zeros(length(nID),nMori);
            for k=1:nMori
                %find relationships of type k
                [row]=find(eFType==k);
                eFFamilyIDsub=eFFamilyID(row,:);
                egroupFEffSFsub=eFEffSF(row,:);
                [famList]=unique(eFFamilyIDsub);
                for j=1:length(famList)
                    val=mean(egroupFEffSFsub(eFFamilyIDsub==famList(j)));
                    nCEffSF(nCFamilyID==famList(j),k)=val;
                    nEffSF(nFamilyID==famList(j),k)=val;
                end
                G_clust.Nodes.EffSF(nCID,k)=nCEffSF(:,k);
                G_frag.Nodes.EffSF(nID,k)=nCEffSF(:,k);
            end 
        end
    end    
end