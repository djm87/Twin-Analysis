function [G_Complete,runCleanupAgain ]= CleanFamilyTree(G_Complete,grains)
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

    G_Complete.Nodes.isTwin = zeros(length(G_Complete.Nodes.Id),1);
    G_Complete.Nodes.isAParent = zeros(length(G_Complete.Nodes.Id),1,'logical');
    G_Complete.Edges.Parent = zeros(size(G_Complete.Edges.pairs),'logical');

    runCleanupAgain=false;
    %loop over groups
    for i=1:max(G_Complete.Edges.Group) 
        egroupId = find((i==G_Complete.Edges.Group)==true); %converts logical arrays to indices
        ngroupId = find((i==G_Complete.Nodes.Group)==true);
        nFamily = G_Complete.Nodes.FamilyID(ngroupId);
        nId = G_Complete.Nodes.Id(ngroupId);
        eType = G_Complete.Edges.type(egroupId);
        eVote = G_Complete.Edges.Vote(egroupId,:);
        ePairs = G_Complete.Edges.pairs(egroupId,:);
        eFamily = G_Complete.Edges.FamilyID(egroupId,:);
        eGlobalId = G_Complete.Edges.GlobalID(egroupId);
        %Initialize parent
        Parent = zeros(size(ePairs,1),2,'logical');

        %Make list of each family relation  
        FamilyRelationList = cell(max(nFamily),max(eType));
        for j = 1:max(nFamily)
            for k = 1:max(eType)
            FamilyInPair = ismember(eFamily(:,:),j);
            [r,c]=find(FamilyInPair);
            FamilyInPair(FamilyInPair) = eType(r)==k;
            FamilyRelationList{j,k} = FamilyInPair;
            end
        end
        
        %Check if a family has the same relationship with more than one
        %twin varient
        if max(nFamily)>100%usualy 2 %i.e. more than parent and one twin type
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
        %Sum the votes for a family 
        if any(sum(Parent,2)==2)
            Parent = zeros(size(ePairs,1),2,'logical');
        end
        %To Do:
        %start building the tree here!
        %Add circular 
        %eVotesSummed are the votes for each for a family and a type
        eVotesSummed = zeros(max(nFamily),max(eType));
        for j = 1:max(nFamily)
            for k = 1:max(eType)
                eVotesSummed(j,k) = sum(eVote(FamilyRelationList{j,k}(:,:)));
            end
        end
        [~,parentID] = sort(eVotesSummed,'descend'); %sorts each column
        for j = 1:size(parentID,1)
            for k = 1:max(eType)
                toSet=FamilyRelationList{parentID(j,k),k};
                notSet=sum(Parent,2)~=0;
                toSet([notSet notSet])=false;
                Parent(toSet) = true;
            end
            if any(sum(Parent,2)==2)
                id=find(sum(Parent,2)==2)
                Parent(id,1)=false;
            end
            if all(sum(Parent,2)) 
                break;
            end
        end
        
%         ismember(eFamily(:,:),j);
        
        FamilyMatrix=zeros(max(nFamily),'logical');
        FamilyTree=zeros(max(nFamily),1)
        %Make Family Tree 
        for j = 1:size(eFamily,1)
            p=eFamily(j,Parent(j,:));
            c=eFamily(j,~Parent(j,:)); 
            FamilyMatrix(p,c)=true;
        end
        FamilyMatrix ;
        %Find child of none
        FamilyTreeParent=find(sum(FamilyMatrix,1)==0);
%         assert(~isempty(FamilyTreeParent),'no clear parent from relationships!')
        if 1==0
            %visualize grain to debug             
            figure; 
            plot(grains(nId),...
                G_Complete.Nodes.FamilyID(nId),'Micronbar','off')
%                 grains(nId).meanOrientation,'Micronbar','off')
            hold on
            e2keep=(i==G_Complete.Edges.Group)==true;
            
%             Ggrain=rmedge(G_Complete,G_Complete.Edges.pairs(~e2keep,1),G_Complete.Edges.pairs(~e2keep,2));
            p=plot(G_Complete,'XData',G_Complete.Nodes.centroids(:,1),...
                'YData',G_Complete.Nodes.centroids(:,2),'displayName','graph');
            hold off
            p.Marker='s';p.NodeColor='k';p.MarkerSize=3;p.EdgeColor='k';
            labeledge(p,G_Complete.Edges.pairs(:,1),G_Complete.Edges.pairs(:,2),G_Complete.Edges.GlobalID);
       
        end
        
        
%         assert(any(sum(FamilyMatrix,1)<3),'Circular relationship between more than two grains')
        circularFamily=find(sum(FamilyMatrix,1)>1)
        rEdge=[];    
        if isempty(FamilyTreeParent)
           %circular relationship exists but it is in the parents column!
           circularFamily=1:max(nFamily);
        end
        for j=1:length(circularFamily)
            cF=circularFamily(j);
            pF=find(FamilyMatrix(:,cF));

            %It is possible for circular relationship to exist for a
            %fragment of a Family while for other fragments, it does not exist. 
            %Since we only want to break one edge not the edges for the whole 
            %family, we need to identify the particular edges in the circular relationship
            
            %Family information needs to be combined with the edge
            %information!

                %Start working with edges instead of families
                EdgeMatrix=zeros(length(nId),length(nId),5,'uint8');
                for k = 1:size(ePairs,1)
                    p=find(ePairs(k,Parent(k,:))==nId);
                    c=find(ePairs(k,~Parent(k,:))==nId);
                    EdgeMatrix(p,c,1)=1;
                    EdgeMatrix(p,c,2)=k;
                    try
                    EdgeMatrix(p,c,3)=eFamily(k,Parent(k,:));
                    catch
                       tic 
                    end
                    EdgeMatrix(p,c,4)=eFamily(k,~Parent(k,:));
                end
                for k = 1:length(pF)
                EdgeMatrix(:,:,5)=EdgeMatrix(:,:,5)+uint8(and(EdgeMatrix(:,:,3)==pF(k),EdgeMatrix(:,:,4)==cF));
                end
                circularEdge=find(sum(EdgeMatrix(:,:,5),1)>1);
    %             assert(all(sum(EdgeMatrix,1)<3),'Need to handle multiple circular relationships with single Parent')
            for k=1:length(circularEdge)
                    %In the case of circular twin relations a child has more than one parent. The script 
                    %compares the boundary ratio between the two parents. If one parent has 90%
                    %more boundary, then it is the parent. Otherwise the parent with the
                    %largest absolute schmid value is chosen as the parent 

        %             G.Edges.FRgB(egroupId,:)

    %                 circularEdgeId=


                    cE=circularEdge(k);
                    pE=find(EdgeMatrix(:,cE,1)==1);
                    eId=zeros(length(pE),1);
                    eFId=zeros(length(pE),1);
                    FRgB=zeros(length(pE),1);
                    EffSF=zeros(length(pE),1);
                                      

                    for kk=1:length(pE)
                        eId(kk)=EdgeMatrix(pE(kk),cE,2);
                        eFId(kk)=eFamily(eId(kk),Parent(eId(kk),:));
                        FRgB(kk)=G_Complete.Edges.FRgB(egroupId(eId(kk)),Parent(eId(kk),:));
                        EffSF(kk)=G_Complete.Edges.EffSF(egroupId(eId(kk)),Parent(eId(kk),:));    
                    end
                     
                    eGlobalId(eId);
                    [ePairs(eId,:),eFamily(eId,:),Parent(eId,:)];
                    
                    
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
                            eGlobalId(eId(EffSF~=EffSF_sorted(end)));
                        end
                    end


            end           
        end
        
        eGlobalId(rEdge);
        %Remove the edges
        rEdgeId=egroupId(rEdge);
        removeEdges=zeros(size(G_Complete.Edges.pairs,1),1,'logical');
        removeEdges(rEdgeId)=true;
        G_Complete=rmedge(G_Complete,G_Complete.Edges.pairs(removeEdges,1),...
            G_Complete.Edges.pairs(removeEdges,2));
        %Reinitialize group quantities
        egroupId = find((i==G_Complete.Edges.Group)==true); %converts logical arrays to indices
        eType = G_Complete.Edges.type(egroupId);
        eVote = G_Complete.Edges.Vote(egroupId,:);
        ePairs = G_Complete.Edges.pairs(egroupId,:);
        eFamily = G_Complete.Edges.FamilyID(egroupId,:);
        eGlobalId = G_Complete.Edges.GlobalID(egroupId);
        Parent(rEdge,:)=[];   
        G_Complete.Edges.Parent(egroupId,:)=Parent;

        %Recalculate Parent 
        FamilyMatrix=zeros(max(nFamily),'logical');
        FamilyTree=zeros(max(nFamily),1);
        %Make Family Tree 
        for j = 1:size(eFamily,1)
            p=eFamily(j,Parent(j,:));
            c=eFamily(j,~Parent(j,:)); 
            FamilyMatrix(p,c)=true;
        end
        [ePairs,eFamily,Parent];

        if (~isempty(find(sum(FamilyMatrix,1)>1)));
            runCleanupAgain=true;
        end
%         assert(isempty(find(sum(FamilyMatrix,1)>1)),'the circular relation routine failed this case.. debug')
        
        
%         if length(FamilyTreeParent)>1
%            %We need to split the parent
%            
%         end
        

            
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
%         %% Important 3/5/19
%         TwinLabel=-1.*ones(size(ePairs));
%         TwinLabel(Parent)=0;
%         for j=1:size(ePairs,1)
%             if all(Parent(j,:)==true)
%                 error('edge with parent definition')
%             elseif Parent(j,2)==true
%                 TwinLabel(j,1)=eType(j); 
%             elseif Parent(j,1)==true
%                 TwinLabel(j,2)=eType(j); 
%             end
%         end
%         %% Important 3/5/19
        
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
%         %% Important 3/5/19
%         for j=1:length(nId)
%             loc1=find(ePairs(:,1)==nId(j));
%             loc2=find(ePairs(:,2)==nId(j));
%             
%             %find a row that has the parent
%             if ~isempty(loc1)
%                 [row,col] = find(TwinLabel(loc1,:)==0,1);
%                 if ~isempty(row)
%                     G_Complete.Nodes.Type(ngroupId(j))=TwinLabel(loc1(row),1);
%                 else
%                     G_Complete.Nodes.Type(ngroupId(j))=-1;
%                 end
%             end
%             if ~isempty(loc2)
%                 [row,col] = find(TwinLabel(loc2,:)==0,1);
%                 if ~isempty(row)
%                     G_Complete.Nodes.Type(ngroupId(j))=TwinLabel(loc2(row),2);
%                 else
%                     G_Complete.Nodes.Type(ngroupId(j))=-1;
%                 end
%             end
%         end
%         %% Important 3/5/19

    end
end
