function [G_Clean,time]= CleanFamilyTree(G,grains,mergedGrains,twin,seg_angle_grouped,runCleanup,maxIter,time)
    
    tic

    CleanupIter=0;
    time.EdgeParents=0;
    G_Clean=G;
    while runCleanup && CleanupIter<maxIter
        tic
        [G_Clean,runCleanup] = EdgeParents(G_Clean,grains,mergedGrains,twin);
        time.EdgeParents=time.EdgeParents+toc;
        if runCleanup
            %Evaluate new clusters in case a merged grain gets seperated
            [G_Clean,time]= AssignFamilyIDs(G_Clean,grains,mergedGrains,seg_angle_grouped,twin,false,time);

        end

        CleanupIter=CleanupIter+1;
    end
    if CleanupIter == maxIter
        disp('More CleanupIter than expected, likely need to solve something manually or improve the code')
    end
    
    if ~isfield(time,'CleanFamilyTree')
        time.CleanFamilyTree=0;
    end
    time.CleanFamilyTree=time.CleanFamilyTree+toc;
end

function [G,runCleanupAgain]= EdgeParents(G,grains,mergedGrains,twin)
%MakeFamilyTree the base parent and all twins that stem from the parent
%
%The parent (1st gen)- a child of no other families
%The parent (2+ gen)- a child of the parent (1st gen) and parent of a
%higher order twin

%In the case of circular twin relations a child has two parents. The script 
%compares the boundary ratio between the two parents. If one parent has 90%
%more boundary, then it is the parent. Otherwise the parent with the
%largest absolute schmid value is chosen as the parent 
%==========================================================================

    G.Nodes.isTwin = zeros(length(G.Nodes.Id),1);
    G.Nodes.isAParent = zeros(length(G.Nodes.Id),1,'logical');
    G.Edges.Parent = zeros(size(G.Edges.pairs),'logical');

    runCleanupAgain=false;
    openType=length(twin);
    
    %Make sure no unknown grains
    %loop over groups
    groups=unique(G.Edges.Group);
    i=1;
    while i<length(groups)+1
        group=groups(i);
        %load edge and node properties for clustered fragments
        egroupId = find((group==G.Edges.Group)==true); %converts logical arrays to indices
        ngroupId = find((group==G.Nodes.Group)==true);
%         typeKnown=~G.Nodes.typeUnknown(ngroupId);
        nId = G.Nodes.Id(ngroupId);
        nFamily = G.Nodes.FamilyID(ngroupId);
        eType = G.Edges.type(egroupId);
        eVote = G.Edges.Vote(egroupId,:);
        ePairs = G.Edges.pairs(egroupId,:);
        eFamily = G.Edges.FamilyID(egroupId,:);
        eGlobalId = G.Edges.GlobalID(egroupId);        
        if ~isempty(eType)
            
            %Filter out nodes and edges that only have type 5 associated 
            %with them. This is important since type 5 don't get a family
            %generation associated with them
            nId_known=unique(ePairs(eType~=openType,:));
            nId_mixed=unique(ePairs(eType==openType,:));
            for j=1:length(nId_mixed)
                if ~any(nId_mixed(j)==nId_known)
                    nInd=find(nId_mixed(j)==nId);
                    nId(nInd)=[];
                    nFamily(nInd)=[];
                end
            end
            etypeKnown=eType~=openType;
            egroupId=egroupId(etypeKnown);
            eType = G.Edges.type(egroupId);
            eVote = G.Edges.Vote(egroupId,:);
            ePairs = G.Edges.pairs(egroupId,:);
            eFamily = G.Edges.FamilyID(egroupId,:);
            eGlobalId = G.Edges.GlobalID(egroupId);  
            
            %renumber families that are around 
            list=unique(nFamily);
            for j=1:length(list)
               nFamily(nFamily==list(j))=j;
               eFamily(eFamily==list(j))=j;
            end
            
            %Returns max(nFamily) x max(eType) cell array containing a logical
            %array of size(ePairs) of the family correlated with the edge Type 
            FamilyRelationList = FamilyRelationships(nFamily,eType,eFamily); 

            %eVotesSummed are the votes for each for a family and a type
            [eVotesSummed,parentID] = buildeVotesSummed(nFamily,eType,eVote,FamilyRelationList);

            %Determine the parent by vote
            Parent = parentByVote(FamilyRelationList,parentID,eType,ePairs,eGlobalId);   
            G.Edges.Parent(egroupId,:)=Parent;
            %Make Family relation matrix 
       
            FamilyMatrix = buildFamilyMatrix(nFamily,eFamily,Parent);

            %Find child of none (i.e. the parent)
            ParentOfAll=find(sum(FamilyMatrix,1)==0);
            
            %Find child with potentially more than one parent
            CircularRelationshipExists=any(sum(FamilyMatrix,1)>1);
            
            %Handle the various problem cases
            if length(ParentOfAll) > 1
            fprintf('More than one parent is apparent.. fix manually\n')
            FamilyMatrix
            [G,runCleanupAgain,i,exitCleanFamily] = ClusterEditor(group,G,grains,mergedGrains,grains.meanOrientation,i,false,false);
            if exitCleanFamily
                break;
            end
            
            elseif isempty(ParentOfAll) || CircularRelationshipExists
               [rEdge,rEdgeGlobalId] = fixCircularFamily(G,FamilyMatrix,egroupId,ePairs,eFamily,nId,nFamily,Parent,eGlobalId);
                %Remove the edges
                G=removeEdge(G,rEdge,egroupId);
                
                %Store remove edges in case we want to look at them later 
                if ~isempty(rEdge)
                    fid = fopen('eRemoveList.txt', 'a+');
                    for j=1:length(rEdge)
                        fprintf(fid, '%d\n', eGlobalId(rEdge(j)));
                    end
                    fclose(fid);
                end
                
                %Reinitialize group quantities
                egroupId = find((group==G.Edges.Group)==true); %converts logical arrays to indices
                eType = G.Edges.type(egroupId);
                etypeKnown=eType~=openType;
                eType=eType(etypeKnown);
                egroupId=egroupId(etypeKnown);
                eVote = G.Edges.Vote(egroupId,:);
                ePairs = G.Edges.pairs(egroupId,:);
                eFamily = G.Edges.FamilyID(egroupId,:);
                eGlobalId = G.Edges.GlobalID(egroupId);   
                Parent(rEdge,:)=[];   
                G.Edges.Parent(egroupId,:)=Parent;

                %Remake Family matrix 
                FamilyMatrix = buildFamilyMatrix(nFamily,eFamily,Parent);

                %Determine if we need to run the cleanup another time
                if ~isempty(find(sum(FamilyMatrix,1)>1))
%                     if isempty(rEdge)
%                         fprintf('There may still be a circular relationship.. inspect manually\n')
%                         FamilyMatrix
%                         [G,runCleanupAgain,~] = ClusterEditor(group,G,grains,mergedGrains,grains.meanOrientation,i,false,false)
%                     else
%                         runCleanupAgain=true;
%                     end
                end
                
                i=i+1;
            else
                G.Edges.Parent(egroupId,:)=Parent;
                i=i+1;
            end
        else
            i=i+1;
        end
    end
end
function G=removeEdge(G,rEdge,egroupId)
    rEdgeId=egroupId(rEdge);
    removeEdges=zeros(size(G.Edges.pairs,1),1,'logical');
    removeEdges(rEdgeId)=true;
    G=rmedge(G,G.Edges.pairs(removeEdges,1),...
        G.Edges.pairs(removeEdges,2));
end
function EdgeMatrix = buildEdgeMatrix(pF,cF,nId,ePairs,eFamily,Parent)
    EdgeMatrix=zeros(length(nId),length(nId),5,'uint8');
    for k = 1:size(ePairs,1)
        p=find(ePairs(k,Parent(k,:))==nId);
        c=find(ePairs(k,~Parent(k,:))==nId);
        EdgeMatrix(p,c,1)=1;
        EdgeMatrix(p,c,2)=k;
%         try
        EdgeMatrix(p,c,3)=eFamily(k,Parent(k,:));
%         catch
%            tic 
%         end
        EdgeMatrix(p,c,4)=eFamily(k,~Parent(k,:));
    end
    for k = 1:length(pF)
        EdgeMatrix(:,:,5)=EdgeMatrix(:,:,5)+uint8(and(EdgeMatrix(:,:,3)==pF(k),EdgeMatrix(:,:,4)==cF));
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
function Parent = parentByVote(FamilyRelationList,parentID,eType,ePairs,eGlobalId)
    %Initialize parent
    Parent = zeros(size(ePairs,1),2,'logical');
    notTwin=load('notTwin.txt');
    notParent=load('notParent.txt');
    eRelation=load('eRelation.txt');
    
    for j=1:length(notTwin)
        Parent(notTwin(j)==ePairs)=true;
    end
    for j=1:length(notParent)
        Parent(fliplr(notParent(j)==ePairs))=true;
    end
    for j=1:size(eRelation,1)
        eId=eRelation(j,1)==eGlobalId;
        if any(eId)
            Parent([eId,eId])=eRelation(j,2:3);
        end
    end
    
    for j = 1:size(parentID,1)
            for k = 1:max(eType)
                toSet=FamilyRelationList{parentID(j,k),k};
                notSet=sum(Parent,2)~=0;
                toSet([notSet notSet])=false;
                Parent(toSet) = true;
            end
            if any(sum(Parent,2)==2)
                id=find(sum(Parent,2)==2) %How do we avoid this altogether?
                Parent(id,1)=false;
            end
            if all(sum(Parent,2)) 
                break;
            end
    end
    
end
function FamilyMatrix = buildFamilyMatrix(nFamily,eFamily, Parent)
    FamilyMatrix=zeros(max(nFamily),'logical');
    %Make Family Tree 
    for j = 1:size(eFamily,1)
        p=eFamily(j,Parent(j,:));
        c=eFamily(j,~Parent(j,:)); 
        FamilyMatrix(p,c)=true;
    end
end
function [rEdge,rEdgeGlobalId] = fixCircularFamily(G,FamilyMatrix,egroupId,ePairs,eFamily,nId,nFamily,Parent,eGlobalId)
    %circular relationship exists but it is in the parents column!
    rEdge=[]; rEdgeGlobalId=[];
    
    circularFamily=1:max(nFamily);
    for j=1:length(circularFamily)
        cF=circularFamily(j);
        pF=find(FamilyMatrix(:,cF));

        %It is possible for circular relationship to exist for a
        %fragment of a Family while for other fragments, it does not exist. 
        %Since we only want to break one edge not the edges for the whole 
        %family, we need to identify the particular edges in the circular relationship

        %Family information needs to be combined with the edge
        %information!
        EdgeMatrix = buildEdgeMatrix(pF,cF,nId,ePairs,eFamily,Parent);
        circularEdge=find(sum(EdgeMatrix(:,:,5),1)>1);
        
        for k=1:length(circularEdge)
            %In the case of circular twin relations a child has more than one parent. The script 
            %compares the boundary ratio between the two parents. If one parent has 90%
            %more boundary, then it is the parent. Otherwise the parent with the
            %largest absolute schmid value is chosen as the parent 

            cE=circularEdge(k);
            pE=find(EdgeMatrix(:,cE,1)==1);
            eId=zeros(length(pE),1);
            eFId=zeros(length(pE),1);
            FRgB=zeros(length(pE),1);
            EffSF=zeros(length(pE),1);


            for kk=1:length(pE)
                eId(kk)=EdgeMatrix(pE(kk),cE,2);
                eFId(kk)=eFamily(eId(kk),Parent(eId(kk),:));
                FRgB(kk)=G.Edges.FRgB(egroupId(eId(kk)),Parent(eId(kk),:));
                EffSF(kk)=G.Edges.EffSF(egroupId(eId(kk)),Parent(eId(kk),:));    
            end


            %group edges with same family edge relationship and
            %assign max schmid for the family
            [uniqueParent,IA,IC]=unique(eFId);
            vecUnique=1:length(uniqueParent);
            if length(uniqueParent)~=1
                for kk=vecUnique
                   group=kk==IC;
                   EffSF(group)=max(EffSF(group));
                end

                %Vote on boundary
                [FRgBls,IA,IC]=unique(FRgB);
                pdiff=abs(FRgBls(end)-FRgBls(1))/((FRgBls(end)+FRgBls(1))/2);

                %Vote on the parent
                %if the same relationship as the last edge, vote
                %the same way
                if pdiff>0.9 
                    rEdge=[rEdge;eId(FRgB~=FRgBls(end))];
                    eGlobalId(eId(FRgB~=FRgBls(end)));

                else
                    [EffSF_sorted,EffSF_I]=sort(EffSF);
                    rEdge=[rEdge;eId(EffSF~=EffSF_sorted(end))];
                    rEdgeGlobalId=eGlobalId(eId(EffSF~=EffSF_sorted(end)));
                end
            end
        end           
    end
end
