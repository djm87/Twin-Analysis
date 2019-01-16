G_Complete.Nodes.isTwin=zeros(length(G_Complete.Nodes.Id),1);
G_Complete.Nodes.isAParent=zeros(length(G_Complete.Nodes.Id),1,'logical');
for i=1:max(G_Complete.Edges.Group)
    egroupId= find((i==G_Complete.Edges.Group)==true); %converts logical arrays to indices
    ngroupId= find((i==G_Complete.Nodes.Group)==true);
    nFamily = G_Complete.Nodes.FamilyID(ngroupId)
    nId = G_Complete.Nodes.Id(ngroupId)
    eType = G_Complete.Edges.type(egroupId)
    eVote = G_Complete.Edges.Vote(egroupId,:)
    ePairs = G_Complete.Edges.pairs(egroupId,:)
    eFamily = G_Complete.Edges.FamilyID(egroupId,:)
    
    %Initialize parent
    Parent = zeros(size(ePairs,1),2,'logical');
    
    %Make list of each family relation  
    FamilyRelationList=cell(max(nFamily),1);
    for j=1:max(nFamily)
        FamilyInPair=ismember(eFamily(:,:),j)
        FamilyRelationList{j}=FamilyInPair
        FamilyRelationList{j}
        eFamily
    end
    
    %The number of relation types for each family 
    %Here we also assign the parent for the case of two tensile variants 
    %having the same family.
    numTypeFamily=zeros(max(nFamily),1);
    for j=1:max(nFamily)
        id=logical(sum(ismember(eFamily(:,:),j),2))
        uniqueTypes=unique(eType(id))
        numTypeFamily(j)=length(uniqueTypes)
        if numTypeFamily(j)>1
            if ismember([1,2],uniqueTypes) %tensile twins can't make tensile twins
                Parent(FamilyRelationList{j})=true;
            end
        end
    end
    
    %Assign the rest of the parents
    for j=1:size(Parent,1)
        if sum(Parent(j,:))==0
            [~,loc]=max(eVote(j,:));
            Parent(j,loc)=true;
        end
    end
    
    %Get twin label for each twin parent
    %Need to take care of circular relationships
    nType=zeros(length(nId),1);
    for j=1:length(nId)
        %determine the grain Types
        %twin nodes 
        if any(nId(j)==ePairs(~Parent))
            nType(j)=unique(eType(nId(j)==ePairs(~Parent)))
        else
            nType(j)=0;
        end
    end
    
    %Finally assign twin label to nodes (need to differentiate generation)
    G_Complete.Nodes.isTwin(ngroupId(ismember(nId,ePairs(~Parent))))=nType(ismember(nId,ePairs(~Parent))); 
    G_Complete.Nodes.isAParent(ngroupId(ismember(nId,ePairs(Parent))))=true ;
end

figure; 
plot(grains(G_Complete.Nodes.Id),G_Complete.Nodes.isTwin,'Micronbar','off')
figure; 
plot(grains(G_Complete.Nodes.Id),grains(G_Complete.Nodes.Id).meanOrientation,'Micronbar','off')


