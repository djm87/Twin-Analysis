function [G_Family,G_clust,G] = FamilyGraph(G_clust,G,grains,computemGrainId,opt)
    %FamilyGraph creates a graph for each merged grain in ClusterGraph
    %The considers the various combinations of relationships between
    %familes

    %Get groups to compute from G_clust 
    toCompute=G_clust.Nodes.isNewGroup | G_clust.Nodes.computeFamily;
    ind=intersect_wrepeat(computemGrainId,G_clust.Nodes.Group);
    toCompute(ind)=true;
    groupList=unique(G_clust.Nodes.Group(toCompute));
    groups=G_clust.Nodes.Group;
    nGroups=length(groupList);
    
    if isempty(groupList)
       fprintf('Grouplist empty, No computation to be done\n')
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
    [G_Family,nEdges] = InitializeFamilyGraph(pairSize,pairIndexes,nFamilyGroup,groupList,nGroups);
    
    % Construct Family twin relationship graph 
    [oriFamily,mori,eFType] = getOriFamily(G_Family,G_clust,grains,groupList,nGroups,nEdges,opt);

    % Test Family pairs for twin relationship
    tol=zeros(opt.nTwin,1);
    tolRlx=zeros(opt.nTwin,1);
    for i=1:opt.nTwin
        tol(i)=opt.twin{i}.tol.misMean;
        tolRlx(i)=opt.twin{i}.tol.misMeanRlx;
    end
    [~,G_Family.Edges.meanType] = TestTwinRelationship(mori,tol,opt,eFType);
    [~,G_Family.Edges.meanTypeRlx] = TestTwinRelationship(mori,tolRlx,opt,eFType); 
    
    %Compute Schmid info for twin/parents in families
    [G_Family,edgeList] = GetSchmidRelative(G_Family,groupList,oriFamily,G_Family.Edges.meanTypeRlx,opt);
    
    %Perform family votes
    [G_Family] = FamilyVotes(G_Family,G_clust,groupList,grains,opt);    

    %Get the parent relationships based on votes
    [G_Family] = getParent(groupList,G_Family) 
    
    %Store the schmid factor in graphs G and G_clust for plotting
    [G_clust,G] = StoreSchmid(G_Family,G_clust,G,groupList,opt);
    
    
    %Update isNewGroup so computation doesn't happen again.
    G_clust.Nodes.isNewGroup(:)=false;
    G.Nodes.isNewGroup=G_clust.Nodes.isNewGroup;

end
function [G_Family] = getParent(groupList,G_Family) 
    G_Family.Edges.Parent = zeros(size(G_Family.Edges.pairs),'logical');    
    
    for i=1:length(groupList)
        group=groupList(i);
        %load edge and node properties for Family graph
        G_Family_sub = subgraph(G_Family,find(group==G_Family.Nodes.Group));
        nFamily = G_Family_sub.Nodes.Family;
        eType = G_Family_sub.Edges.meanTypeRlx;
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
        G_Family.Edges.Parent(eGlobalID,:)=Parent; %Need to add back in
    end
    %update the directions of the graph and remove non-twin edges
    toRemove=find(all(G_Family.Edges.Parent==0,2));
    G_Family=rmedge(G_Family,toRemove);
    G_Family.Edges.eGlobalID=[1:numedges(G_Family)]';
    toFlip=G_Family.Edges.Parent(:,2);
    G_Family.Edges.pairs(toFlip,:)=fliplr(G_Family.Edges.pairs(toFlip,:));
    G_Family.Edges.FamilyID(toFlip,:)=fliplr(G_Family.Edges.FamilyID(toFlip,:));
    G_Family.Edges.eRelationship(toFlip,:)=fliplr(G_Family.Edges.eRelationship(toFlip,:));
    G_Family.Edges.sigma13(toFlip,:)=fliplr(G_Family.Edges.sigma13(toFlip,:));
    G_Family.Edges.FRArea(toFlip,:)=fliplr(G_Family.Edges.FRArea(toFlip,:));
    G_Family.Edges.FRgB(toFlip,:)=fliplr(G_Family.Edges.FRgB(toFlip,:));
    G_Family.Edges.FREffSF(toFlip,:)=fliplr(G_Family.Edges.FREffSF(toFlip,:));
    G_Family.Edges.Vote(toFlip,:)=fliplr(G_Family.Edges.Vote(toFlip,:));
    G_Family.Edges.Parent(toFlip,:)=fliplr(G_Family.Edges.Parent(toFlip,:));
    G_Family=flipedge(G_Family,find(toFlip));
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
function [G_Family,nEdges] = InitializeFamilyGraph(pairSize,pairIndexes,nFamilyGroup,groupList,nGroups)

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
    G_Family=digraph(s,t);
    G_Family.Edges.pairs=[s,t];
    G_Family.Edges.Group=Group;
    G_Family.Edges.FamilyID=FamilyID;
    G_Family.Edges.eRemove=zeros(nEdges,1,'logical');
    G_Family.Edges.eIsParentAll=zeros(nEdges,1,'logical');
    G_Family.Edges.eNotParentAll=zeros(nEdges,1,'logical');
    G_Family.Edges.eNotParent=zeros(nEdges,1,'logical');
    G_Family.Edges.eRelationship=zeros(nEdges,2,'int8');
    G_Family.Edges.eGlobalID=[1:nEdges]';
    G_Family.Nodes.GlobalID=[1:numnodes(G_Family)]';
    % Construct node Family and Group
    [~,ind]=unique(G_Family.Edges.pairs);
    G_Family.Nodes.Family=FamilyID(ind);
    nGroup(s)=Group;
    nGroup(t)=Group;
    G_Family.Nodes.Group=nGroup';
end
function [oriFamily,mori,eFType] = getOriFamily(G_Family,G_clust,grains,groupList,nGroups,nEdges,opt)
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
        egroupFId = find(group==G_Family.Edges.Group);
        nFamily = G_clust.Nodes.FamilyID(ngroupId);     
        eFamily = G_clust.Edges.FamilyID(egroupId,:);
        eType = G_clust.Edges.type(egroupId);
        FamilyID=G_Family.Edges.FamilyID(egroupFId,:);

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
    
    mori=inv(ori1).*ori2; 
    oriFamily=[ori1,ori2];
end
function [G_clust,G] = StoreSchmid(G_Family,G_clust,G,groupList,opt)
    
    %First store new info in G_clust - EffSF
    nTwin=opt.nTwin;
    for i=1:length(groupList) 
        group=groupList(i);
        eFID = find(group==G_Family.Edges.Group);
        if ~isempty(eFID)
%             nFID = find(group==G_Family.Nodes.Group);
%             eCID = find(group==G_clust.Edges.Group);
            nCID = find(group==G_clust.Nodes.Group);            
            nID = find(group==G.Nodes.Group);

            eFEffSF = G_Family.Edges.EffSF(eFID,:);
            eFType = G_Family.Edges.meanTypeRlx(eFID);
            eFFamilyID=G_Family.Edges.FamilyID(eFID,:);
            nCFamilyID = G_clust.Nodes.FamilyID(nCID);
            nFamilyID = G.Nodes.FamilyID(nID);
            
            nCEffSF=zeros(length(nCID),nTwin);
            nEffSF=zeros(length(nID),nTwin);
            for k=1:nTwin
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
                G.Nodes.EffSF(nID,k)=nCEffSF(:,k);
            end 
        end
    end    
end