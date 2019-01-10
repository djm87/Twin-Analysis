function G_Complete = MakeFamilyTree(G_Complete)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    G_Complete.Nodes.isTwin=zeros(length(G_Complete.Nodes.Id),1);
    G_Complete.Nodes.isAParent=zeros(length(G_Complete.Nodes.Id),1,'logical');
    for i=1:max(G_Complete.Edges.Group) %loops over groups
        egroupId= find((i==G_Complete.Edges.Group)==true); %converts logical arrays to indices
        ngroupId= find((i==G_Complete.Nodes.Group)==true);
        nFamily = G_Complete.Nodes.FamilyID(ngroupId);
        nId = G_Complete.Nodes.Id(ngroupId);
        eType = G_Complete.Edges.type(egroupId);
        eVote = G_Complete.Edges.Vote(egroupId,:);
        ePairs = G_Complete.Edges.pairs(egroupId,:);
        eFamily = G_Complete.Edges.FamilyID(egroupId,:);

        %Initialize parent
        Parent = zeros(size(ePairs,1),2,'logical');

        %Make list of each family relation  
        FamilyRelationList=cell(max(nFamily),1);
        for j=1:max(nFamily)
            FamilyInPair=ismember(eFamily(:,:),j);
            FamilyRelationList{j}=FamilyInPair;
        end
        
        %Sum the votes for a family
        eVotesSummed=zeros(max(nFamily),1);
        for j=1:max(nFamily)
            eVotesSummed(j)=sum(eVote(FamilyRelationList{j}(:,:)));
        end
        parentID=find(eVotesSummed==max(eVotesSummed));
        Parent(FamilyRelationList{parentID})=true;

        %The number of relation types for each family 
        %Here we also assign the parent for the case of two tensile variants 
        %having the same family.
%         numTypeFamily=zeros(max(nFamily),1);
%         for j=1:max(nFamily)
%             id=logical(sum(ismember(eFamily(:,:),j),2));
%             uniqueTypes=unique(eType(id));
%             numTypeFamily(j)=length(uniqueTypes);
%             if numTypeFamily(j)>1
%                 if ismember([1,2],uniqueTypes) %Tensile Twin can't make tensile twin
%                     Parent(FamilyRelationList{j})=true;
%                 end
%             end
%         end

        %Assign pairwise labels
        TwinLabel=-1.*ones(size(ePairs));
        TwinLabel(Parent)=0;
        for j=1:size(ePairs,1)
            if all(Parent(j,:)==true)
                error('edge with parent definition')
            elseif Parent(j,2)==true
                TwinLabel(j,1)=eType(j); 
            elseif Parent(j,1)==true
                TwinLabel(j,2)=eType(j); 
            end
        end
        
        %Here take care of second order twins 
%         loopCnt=0;
%         while any(TwinLabel==0)
%            loopCnt=loopCnt+1;
%             if loopCnt>5
%                 break;
%             end
%             
%             for j=1:size(TwinLabel,1)
%                 if TwinLabel(j,1)==0
%                    tmp=TwinLabel(find(ePairs(:,1)==ePairs(j,1)))
%                    if tmp>0 
%                        TwinLabel(j,1)=max(tmp)
%                    else
%                        tmp=TwinLabel(find(ePairs(:,2)==ePairs(j,1)))
%                        if tmp>0 
%                           TwinLabel(j,1)=max(tmp)
%                        end
%                    end
%                 end
%                 if TwinLabel(j,2)==0
%                    tmp=TwinLabel(find(ePairs(:,2)==ePairs(j,2)))
%                    if max(tmp)>0 
%                        TwinLabel(j,2)=max(tmp)
%                    else
%                        tmp=TwinLabel(find(ePairs(:,1)==ePairs(j,2)))
%                        if max(tmp)>0 
%                           TwinLabel(j,2)=max(tmp)
%                        end
%                    end
%                 end
%             end
%         end
% %         
        %Finally assign twin label to nodes
%         G_Complete.Nodes.isTwin(ngroupId(ismember(nId,ePairs(~Parent))))=5; 

        for j=1:length(nId)
            loc1=find(ePairs(:,1)==nId(j));
            loc2=find(ePairs(:,2)==nId(j));
            
            %find a row that has the parent
            if ~isempty(loc1)
                [row,col] = find(TwinLabel(loc1,:)==0,1);
                if ~isempty(row)
                    G_Complete.Nodes.Type(ngroupId(j))=TwinLabel(loc1(row),1);
                else
                    G_Complete.Nodes.Type(ngroupId(j))=-1;
                end
            end
            if ~isempty(loc2)
                [row,col] = find(TwinLabel(loc2,:)==0,1);
                if ~isempty(row)
                    G_Complete.Nodes.Type(ngroupId(j))=TwinLabel(loc2(row),2);
                else
                    G_Complete.Nodes.Type(ngroupId(j))=-1;
                end
            end
        end
    end
end
