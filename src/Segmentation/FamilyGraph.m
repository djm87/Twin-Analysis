function [G_Family,G_clust] = FamilyGraph(G_clust,grains,mGrains,opt)
    %FamilyGraph creates a graph for each merged grain in ClusterGraph
    %The considers the various combinations of relationships between
    %familes
    
    %Create the minimum set of pairs that should be tested for twins
    %relationships for each family group size
    pairIndexes=cell(max(G_clust.Nodes.FamilyID),1);
    pairSize=zeros(max(G_clust.Nodes.FamilyID),1);
    for numFamily=1:max(G_clust.Nodes.FamilyID)
        [row,col]=find(triu(ones(numFamily))-eye(numFamily));
        pairIndexes{numFamily}=[row,col];
        pairSize(numFamily)=length(row);
    end
    
    % Construct Family twin relationship graph 
    numFamilypGroup=zeros(length(mGrains),1);
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
    G_Family.Edges.eRemove=zeros(length(G_Family.Edges.Group),1,'logical');
    G_Family.Edges.eIsParentAll=zeros(length(G_Family.Edges.Group),1,'logical');
    G_Family.Edges.eNotParentAll=zeros(length(G_Family.Edges.Group),1,'logical');
    G_Family.Edges.eNotParent=zeros(length(G_Family.Edges.Group),1,'logical');
    G_Family.Edges.eRelationship=zeros(length(G_Family.Edges.Group),2,'int8');
    
    % Construct node Family and Group
    [~,ind]=unique(G_Family.Edges.pairs);
    G_Family.Nodes.Family=familyPair(ind);
    nGroup(s)=Group;
    nGroup(t)=Group;
    G_Family.Nodes.Group=nGroup';
    
    % Construct Family twin relationship graph 
    oriAll=grains.meanOrientation ;
    areaAll=grains.area;
    ori1 = orientation.byEuler(zeros(numEdges,1),...
        zeros(numEdges,1),zeros(numEdges,1),...
            'ZYZ',opt.CS{2});
    ori2=ori1;
    eFType=zeros(numEdges,1);
    for i=1:length(mGrains)
        ngroupId = find(i==G_clust.Nodes.Group);
        egroupId = find(i==G_clust.Edges.Group); %converts logical arrays to indices
        egroupFId = find(i==G_Family.Edges.Group);
        nFamily = G_clust.Nodes.FamilyID(ngroupId);     
        eFamily = G_clust.Edges.FamilyID(egroupId,:);
        eType = G_clust.Edges.type(egroupId);
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
    oriFamily=[ori1,ori2];
    tol=zeros(opt.nTwin,1);
    tolRlx=zeros(opt.nTwin,1);
    for i=1:opt.nTwin
        tol(i)=opt.twin{i}.tol.misMean;
        tolRlx(i)=opt.twin{i}.tol.misMeanRlx;
    end
    [~,G_Family.Edges.meanType] = TestTwinRelationship(mori,tol,opt,eFType);
    [~,G_Family.Edges.meanTypeRlx] = TestTwinRelationship(mori,tolRlx,opt,eFType);
    
    %Compute Schmid info for twin/parents in clustered grains
    %This computes Schmid factor for twin/parent identification and 
    %is stored at both edge and node level
    [G_Family,G_clust]=GetSchmidRelative(G_Family,G_clust,oriFamily,G_Family.Edges.meanTypeRlx,grains,mGrains,opt)
    
    %Perform family votes
    [G_Family,G_clust] = FamilyVotes(G_Family,G_clust,unique(G_clust.Nodes.Group),grains,mGrains,opt);    

    %Update isNewGroup so computation doesn't happen again.
    
end
