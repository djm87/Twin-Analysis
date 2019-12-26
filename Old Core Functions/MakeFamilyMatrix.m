function [FamilyMatrix,FamilyRelationList,eVotesSummed]=...
    MakeFamilyRelationship(nId,nFamily,eType,eVote,eFamily,ePairs,useVariantGroup,nVariantGroup)
%MakeFamilyMatrix applies the voting scheme to determine parent and twin
%relationships based on the current clusters of grains. Likely there will
%be some fragments that have more than one parent or there will be
%initially grains groupped that shouldn't be grouped.

%=========================================================================%
%A family in a group is a parent automatically if it has the same twin mode 
%relationship with two or more other familes in the group. This is likely 
%true for monotonic loading only and can be turned off with the flag ADD!!!
%Note that this assertion could be true for a particular twin type or . 
%
 
%==========================================================================
        
        %Initialize the Parent array
        Parent = zeros(size(ePairs,1),2,'logical');

        %Correlate the eType and FamilyID
        %Returns max(nFamily) x max(eType) cell array containing a logical
        %array of size(ePairs) of the family correlated with the edge Type 
        FamilyRelationList = getFamilyRelationList(nFamily,eType,eFamily) 
        
        %Assign parent based on a grain containing multiple twin variants
        if and(max(nFamily)>nVariantGroup,useVariantGroup)
            Parent = parentByVariants(FamilyRelationList,nFamily,eType,eFamily)
        end       

        %eVotesSummed are the votes for each for each family and for each eType
        [eVotesSummed,parentID] = geteVotesSummed(nFamily,eType,eVote,FamilyRelationList);
        
        %Returns a full Parent list, completing the work started by 
        %parentByVariants if it was called
        Parent = parentByVote(FamilyRelationList,Parent,parentID,eType);
        
        %Returns a matrix containing the rows and columns relating parent
        %and child relationships
        FamilyMatrix = buildFamilyMatrix(nFamily,eFamily, Parent);
        EdgeMatrix = buildEdgeMatrix(nId,ePairs,eFamily,Parent)
        
        %Check the parent list isn't messed up 
        assert(all(sum(Parent,2)<2),'both nodes of an edge are the parent')
        
        %Assign parents to global 
        %G_Complete.Edges.Parent(eGlobalId)=Parent
end
function FamilyRelationList = getFamilyRelationList(nFamily,eType,eFamily) 
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
function Parent = parentByVariants(FamilyRelationList,nFamily,eType,eFamily)
    for j = 1:max(nFamily)
        for k = 1:max(eType)
            %Get current edges
            currentEdges=find(sum(FamilyRelationList{j,k},2));
            currentEFamily=eFamily(currentEdges,:);
            currentEType=eType(currentEdges);
            for kk=1:size(currentEFamily,1) 
                for kkk=1:size(currentEFamily,1)
                    if all(currentEFamily(kk,:)==fliplr(currentEFamily(kkk,:)));
                        currentEFamily(kkk,:)=currentEFamily(kk,:);
                    end
                end
            end
            [uniqueFamilies,IA,IC]=unique(currentEFamily,'rows');
            if length(IA)>2
                Parent(FamilyRelationList{j,k}) = true;
            end
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
function [eVotesSummed,parentID] = geteVotesSummed(nFamily,eType,eVote,FamilyRelationList)
    eVotesSummed = zeros(max(nFamily),max(eType));
    for j = 1:max(nFamily)
        for k = 1:max(eType)
            eVotesSummed(j,k) = sum(eVote(FamilyRelationList{j,k}(:,:)));
        end
    end
    [eVotesSummed,parentID] = sort(eVotesSummed,'descend'); %sorts each column
end
function Parent = parentByVote(FamilyRelationList,Parent,parentID,eType)
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
function EdgeMatrix = buildEdgeMatrix(nId,ePairs,eFamily,Parent)
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