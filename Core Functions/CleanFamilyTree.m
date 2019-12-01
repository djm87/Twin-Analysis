function [G,runCleanupAgain ]= CleanFamilyTree(G,grains,tFilter)
%MakeFamilyTree the base parent and all twins that stem from the parent
%
%The parent - a child of no other families
%=========================================================================%
%A family in a group is a parent automatically if it has the same twin mode 
%relationship with two or more other familes in the group. This is likely 
%true for monotonic loading only and can be turned off with the flag ADD!!!
% Note that this assertion could be true for a n generation case. 
%
%In the case of circular twin relations a child has two parents. The script 
%compares the boundary ratio between the two parents. If one parent has 90%
%more boundary, then it is the parent. Otherwise the parent with the
%largest absolute schmid value is chosen as the parent 
%==========================================================================

    G.Nodes.isTwin = zeros(length(G.Nodes.Id),1);
    G.Nodes.isAParent = zeros(length(G.Nodes.Id),1,'logical');
    G.Edges.Parent = zeros(size(G.Edges.pairs),'logical');

    runCleanupAgain=false;
    %loop over groups
    for i=1:max(G.Edges.Group) 
        %load edge and node properties for clustered fragments
        egroupId = find((i==G.Edges.Group)==true); %converts logical arrays to indices
        ngroupId = find((i==G.Nodes.Group)==true);
        nFamily = G.Nodes.FamilyID(ngroupId);
        nId = G.Nodes.Id(ngroupId);
        eType = G.Edges.type(egroupId);
        eVote = G.Edges.Vote(egroupId,:);
        ePairs = G.Edges.pairs(egroupId,:);
        eFamily = G.Edges.FamilyID(egroupId,:);
        eGlobalId = G.Edges.GlobalID(egroupId);        
        
        %Returns max(nFamily) x max(eType) cell array containing a logical
        %array of size(ePairs) of the family correlated with the edge Type 
        FamilyRelationList = FamilyRelationships(nFamily,eType,eFamily); 
             
        %eVotesSummed are the votes for each for a family and a type
        [eVotesSummed,parentID] = buildeVotesSummed(nFamily,eType,eVote,FamilyRelationList);
        
        %Determine the parent by vote
        Parent = parentByVote(FamilyRelationList,parentID,eType,ePairs);   
        
        %Make Family relation matrix 
        FamilyMatrix = buildFamilyMatrix(nFamily,eFamily,Parent);
        
        %Find child of none (i.e. the parent)
        FamilyTreeParent=find(sum(FamilyMatrix,1)==0);
        
        if isempty(FamilyTreeParent)
           [rEdge,rEdgeGlobalId] = fixCircularFamily(G,FamilyMatrix,egroupId,eId,ePairs,eFamily,nId,nFamily,Parent);
            %Remove the edges
            G=removeEdge(G,rEdge,egroupId);
            
            %Reinitialize group quantities
            egroupId = find((i==G.Edges.Group)==true); %converts logical arrays to indices
            eType = G.Edges.type(egroupId);
            eVote = G.Edges.Vote(egroupId,:);
            ePairs = G.Edges.pairs(egroupId,:);
            eFamily = G.Edges.FamilyID(egroupId,:);
            eGlobalId = G.Edges.GlobalID(egroupId);
            Parent(rEdge,:)=[];   
            G.Edges.Parent(egroupId,:)=Parent;
            toc
            %Remake Family matrix 
            FamilyMatrix = buildFamilyMatrix(nFamily,eFamily,Parent);

            %Determine if we need to run the cleanup a second time
            if (~isempty(find(sum(FamilyMatrix,1)>1)));
                runCleanupAgain=true;
            end
        else
            G.Edges.Parent(egroupId,:)=Parent;
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
function Parent = parentByVote(FamilyRelationList,parentID,eType,ePairs)
    %Initialize parent
    Parent = zeros(size(ePairs,1),2,'logical');
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
function [rEdge,rEdgeGlobalId] = fixCircularFamily(G,FamilyMatrix,egroupId,eId,ePairs,eFamily,nId,nFamily,Parent)
    %circular relationship exists but it is in the parents column!
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
        EdgeMatrix = buildEdgeMatrix(pF,cF,nId,ePairs,eFamily,Parent)
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
                    eGlobalId(eId(FRgB~=FRgBls(end)))

                else
                    [EffSF_sorted,EffSF_I]=sort(EffSF);
                    rEdge=[rEdge;eId(EffSF~=EffSF_sorted(end))];
                    rEdgeGlobalId=eGlobalId(eId(EffSF~=EffSF_sorted(end)));
                end
            end
        end           
    end
end
