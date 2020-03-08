function [G_Family,G_clust,exflagGroup]= CleanFamilyTree(groups,G_Family,G_clust,grains,mGrains,exflagGroup,opt)   
%cleanFamilyTree automatically removes circular twin relationships and
%outputs a list of groups that need user input to resolve.

%exflagGroup=0 clean family
%exflagGroup=1 no twin relationships but multiple families
%exflagGroup=2 too many parents
%exflagGroup=3 no parent but twin relationships exist
%exflagGroup=4 max number of generations hit while assigning generation
%exflagGroup=5 not all fragments were related... something is wrong
%exflagGroup=6 single frag
%==========================================================================

    G_Family.Edges.Parent = zeros(size(G_Family.Edges.pairs),'logical');    
   
    %loop over groups
    i=1;
    while i<=length(groups)
        group=groups(i)
        %load edge and node properties for Family graph
        egroupId = find((group==G_Family.Edges.Group)==true); %converts logical arrays to indices
        ngroupId = find((group==G_Family.Nodes.Group)==true);
        nFamily = G_Family.Nodes.Family(ngroupId);
        eType = G_Family.Edges.meanTypeRlx(egroupId);
        eVote = G_Family.Edges.Vote(egroupId,:);
        ePairs = G_Family.Edges.pairs(egroupId,:);
        eFamily = G_Family.Edges.familyPair(egroupId,:);
        nGeneration = -ones(length(ngroupId),1,'int8');
        nType = -ones(length(ngroupId),1,'int8');   
        eRemove = G_Family.Edges.eRemove(egroupId);
        
        %load edge and node properties for cluster graph
        nClustgroupId = find((group==G_clust.Nodes.Group)==true);
        eClustgroupId = find((group==G_clust.Edges.Group)==true); %converts logical arrays to indices
        nClustId = G_clust.Nodes.Id(nClustgroupId);
        nClustFamily = G_clust.Nodes.FamilyID(nClustgroupId);
        eClustFamily = G_clust.Edges.FamilyID(eClustgroupId,:);
        eClustType = G_clust.Edges.type(eClustgroupId);
        nClustGeneration = -ones(length(nClustId),1,'int8');
        nClustType = -ones(length(nClustId),1,'int8');        

        if isempty(eType)
            %Single fragment
            exflagGroup(i)=6;
            i=i+1;
        elseif all(eType==0)
            %Multiple fragments but no relationships
            exflagGroup(i)=1;
            i=i+1;
        else
            %Returns max(nFamily) x max(eType) cell array containing a logical
            %array of size(ePairs) of the family correlated with the edge Type 
            FamilyRelationList = FamilyRelationships(nFamily,eType,eFamily); 

            %eVotesSummed are the votes for each for a family and a type
            [eVotesSummed,parentID] = buildeVotesSummed(nFamily,eType,eVote,FamilyRelationList);

            %Determine the parent by vote
            Parent = parentByVote(eVotesSummed,eType,eFamily);   
            G_Family.Edges.Parent(egroupId,:)=Parent; %Need to add back in
            %Make Family relation matrix 
       
            FamilyMatrix = buildFamilyMatrix(nFamily,eFamily,Parent,eRemove);

            %Find child of none (i.e. the parent) and exclude families with
            %no relationships
            TooManyParentsExist=sum(all(~FamilyMatrix,1) & ~all(~FamilyMatrix,2)')>1;  
            
            %Find child with more than one parent (circular relationship)
            CircularRelationshipExists=any(sum(FamilyMatrix,1)>1);
            
            if TooManyParentsExist
                exflagGroup(i)=2;
                
                if opt.debugFamilyTree
                    fprintf('More than one parent is apparent.. fix manually\n')
                    FamilyMatrix
                    [G_clust,runCleanupAgain,i,exitCleanFamily] = ClusterEditor(group,G_clust,grains,mGrains,grains.meanOrientation,i,1,1,0,1,0); 
                    if exitCleanFamily
                        fprintf('exiting CleanFamilyTree\n')
                        break;
                    end
                else
                   i=i+1;
                end
                
            elseif CircularRelationshipExists
               [rEdge] = fixCircularFamily(G_Family,...
                   FamilyMatrix,egroupId,ePairs,eFamily,eClustFamily,...
                   eClustType,nFamily,Parent,opt);
                for k=1:length(rEdge)
                    pInd=find(Parent(rEdge(k),:));
                    cInd=find(~Parent(rEdge(k),:));
                    FamilyMatrix(eFamily(rEdge(k),pInd),eFamily(rEdge(k),cInd))=false;
                end
                
                %Store remove edges in case we want to look at them later 
                G_Family.Edges.eRemove(egroupId(rEdge))=true;

            else %Everything looks good to process!
                G_Family.Edges.Parent(egroupId,:)=Parent;
                
                %Add Type to the family matrix
                FamilyMatrixType=int8(FamilyMatrix);
                [r,c]=find(FamilyMatrix);
                
                for j=1:length(r)
                    FamilyMatrixType(r(j),c(j))=eType(all(eFamily==[r(j),c(j)],2) | all(fliplr(eFamily)==[r(j),c(j)],2));
                end
                
                %Call labeling function
                exflag=0;
                [nGeneration,nType,exflag] = fillGenerations(...
                    FamilyMatrixType,nGeneration,nType,[],0,8,exflag,opt);
                
                exflagGroup(i)=exflag;
                
                %Store the generation and type data in G_clust for further
                %processing and plotting

                for j=1:max(nFamily)
                    lInd=nClustFamily==j;
                    nClustType(lInd)=nType(j);
                    nClustGeneration(lInd)=nGeneration(j);
                end
                
                G_clust.Nodes.type(nClustgroupId)=nClustType;
                G_clust.Nodes.Generation(nClustgroupId)=nClustGeneration;

                i=i+1;
            end
        end
        G_clust.Nodes.nClustgroupId
    end
end
function G_clust=removeEdge(G_clust,rEdge,egroupId)
    rEdgeId=egroupId(rEdge);
    removeEdges=zeros(size(G_clust.Edges.pairs,1),1,'logical');
    removeEdges(rEdgeId)=true;
    G_clust=rmedge(G_clust,G_clust.Edges.pairs(removeEdges,1),...
        G_clust.Edges.pairs(removeEdges,2));
end
function [nGeneration,nType,exflag] = fillGenerations(FamilyMatrix,...
                                        nGeneration,nType,genId,gen,maxGen,exflag,opt)
    %Simple recursive function to label type and generation 
    %exflag=1 no twin relationships
    %exflag=2 too many parents
    %exflag=3 no parent
    %exflag=4 max number of generations hit
    %exflag=5 not all fragments were related... something is wrong
    if gen==0
        %Check to see if there are twin relationships
        if all(all(FamilyMatrix==0))
            exflag=1;
            return
        end
        
        %See if there is a parent
        genId = find(all(FamilyMatrix==0,1) & ~all(FamilyMatrix==0,2)');
        if numel(genId)==1
            nType(genId)=0;
            nGeneration(genId)=0;
            gen=1;
            
            %Assign unknown relationships a type so they don't get flagged
            %as causing a major failure (i.e see condition for exflag=5)
            nType(all(FamilyMatrix==0,1) &...
                all(FamilyMatrix==0,2)')=opt.twinUnknown;                 
            
        else %No parent or too many parents!
            if numel(genId)>1
                exflag=2;
            else
                exflag=3;
            end
            nType(:)=0;
            nGeneration(:)=-1;
            return
        end
    end
    if gen==maxGen
        %max number of generations hit
        exflag=4;
        
        %Reset nGeneration so we know to debug.
        nGeneration(:)=-1;
        nType(:)=0;
        return
    else
        for i=1:length(genId)
            genId1 = find(FamilyMatrix(genId(i),:)~=0);
            if ~isempty(genId1)
                nGeneration(genId1)=gen;
                nType(genId1)=FamilyMatrix(genId(i),genId1);
                if any(nGeneration==-1)
                    [nGeneration,nType,exflag] = fillGenerations(FamilyMatrix,...
                        nGeneration,nType,genId1,gen+1,maxGen,exflag,opt);
                end
            end
        end
        if any(nGeneration==-1 & nType~=opt.twinUnknown) && gen==1
            if exflag~=4
               exflag=5;
            end
            nGeneration(:)=-1;
            nType(:)=0;
        end
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
    
%     notTwin=load('notTwin.txt');
%     notParent=load('notParent.txt');
%     eRelation=load('eRelation.txt');
%     
%     for j=1:length(notTwin)
%         Parent(notTwin(j)==ePairs)=true;
%     end
%     for j=1:length(notParent)
%         Parent(fliplr(notParent(j)==ePairs))=true;
%     end
%     for j=1:size(eRelation,1)
%         eId=eRelation(j,1)==eGlobalId;
%         if any(eId)
%             Parent([eId,eId])=eRelation(j,2:3);
%         end
%     end
    
%     for j = 1:size(parentID,1)
%             for k = 1:max(eType)
%                 toSet=FamilyRelationList{parentID(j,k),k};
%                 notSet=sum(Parent,2)~=0;
%                 toSet([notSet notSet])=false;
%                 Parent(toSet) = true;
%             end
%             if any(sum(Parent,2)==2)
%                 id=find(sum(Parent,2)==2) %How do we avoid this altogether?
%                 Parent(id,1)=false;
%             end
%             if all(sum(Parent,2)) 
%                 break;
%             end
%     end
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
function FamilyMatrix = buildFamilyMatrix(nFamily,eFamily, Parent,eRemove)
    FamilyMatrix=zeros(max(nFamily),'logical');
    %Make Family Tree 
    for j = 1:size(eFamily,1)
        p=eFamily(j,Parent(j,:));
        c=eFamily(j,~Parent(j,:)); 
        FamilyMatrix(p,c)=~eRemove(j);
    end   
end
function [rEdge] = fixCircularFamily(G_Family,FamilyMatrix,egroupId,ePairs,eFamily,eClustFamily,eClustType,nFamily,Parent,opt)
    
    %Initialize edges
    rEdge=[]; 

    %For each child look for multiple parents
    cFList=find(sum(FamilyMatrix,1)>1);
    
    %Loop over children with too many parents
    for cF=cFList
        %For some child find all the parents
        pF=find(FamilyMatrix(:,cF));

        %In the case of circular twin relations a child has more than 
        %one parent. The script compares the boundary ratio between 
        %the parents. If one parent has 90% more boundary, then it is
        %the parent. Otherwise the parent with the largest absolute 
        %schmid value is chosen as the parent 

        %             cE=circularEdge(k);
        %             pE=find(EdgeMatrix(:,cE,1)==1);
        %Find the right edge for the parent/child relationships
        eId=zeros(length(pF),1);
        FRgB=zeros(length(pF),1);
        EffSF=zeros(length(pF),1);
        localEdges=zeros(length(pF),1);
        
        for kk=1:length(pF)
            PotentialEdges=all(eClustFamily==[pF(kk),cF] | eClustFamily==[cF,pF(kk)],2);
            PotentialEdges(eClustType==opt.twinUnknown)=false;     
            localEdges(kk)=sum(PotentialEdges);
            eId(kk)=find(sum(or(pF(kk)==eFamily,cF==eFamily),2)==2);
            FRgB(kk)=G_Family.Edges.FRgB(egroupId(eId(kk)),Parent(eId(kk),:));
            G_Family.Edges.Vote(egroupId(eId(kk)),Parent(eId(kk),:));
            EffSF(kk)=G_Family.Edges.EffSF(egroupId(eId(kk)),Parent(eId(kk),:));    
        end
        
        if all(localEdges==0)
            %No edges between touching grains. Then use the effective 
            %schmid
                        
            %sort schmid and calculation difference
            [EffSF_sorted,EffSF_I]=sort(EffSF);
            rEdge=[rEdge;eId(EffSF~=EffSF_sorted(end))];
        else
            rEdge=[rEdge;eId(localEdges==0)];
            eId(localEdges==0)=[];
            EffSF(localEdges==0)=[];
            %Check to see if there are still some circular relationships
            if numel(localEdges~=0)>1
                %sort schmid
                [EffSF_sorted,EffSF_I]=sort(EffSF);
                rEdge=[rEdge;eId(EffSF~=EffSF_sorted(end))];
            end
        end
            
    end           
end



            
           
%             elseif CircularRelationshipExists
%                [rEdge,rEdgeGlobalId] = fixCircularFamily(G_Family,...
%                    FamilyMatrix,egroupId,ePairs,eFamily,eClustFamily,...
%                    eClustType,nFamily,Parent,opt);
%                
%                for k=1:length(rEdge)
%                 pInd=find(Parent(rEdge(k),:));
%                 cInd=find(~Parent(rEdge(k),:));
%                 FamilyMatrix(eFamily(rEdge(k),pInd),eFamily(rEdge(k),cInd))=false;
%                end
%                 %Remove the edges
%                 G_clust=removeEdge(G_clust,rEdge,egroupId);
%                 
%                 %Store remove edges in case we want to look at them later 
%                 if ~isempty(rEdge)
%                     fid = fopen('eRemoveList.txt', 'a+');
%                     for j=1:length(rEdge)
%                         fprintf(fid, '%d\n', eGlobalId(rEdge(j)));
%                     end
%                     fclose(fid);
%                 end
%                 
%                 %Reinitialize group quantities
%                 egroupId = find((group==G_clust.Edges.Group)==true); %converts logical arrays to indices
%                 eType = G_clust.Edges.type(egroupId);
%                 etypeKnown=eType~=opt.twinUnknown & eType~=0;
%                 eType=eType(etypeKnown);
%                 egroupId=egroupId(etypeKnown);
%                 eVote = G_clust.Edges.Vote(egroupId,:);
%                 ePairs = G_clust.Edges.pairs(egroupId,:);
%                 eFamily = G_clust.Edges.FamilyID(egroupId,:);
%                 eGlobalId = G_clust.Edges.GlobalID(egroupId);   
%                 Parent(rEdge,:)=[];   
%                 G_clust.Edges.Parent(egroupId,:)=Parent;
% 
%                 %Remake Family matrix 
%                 FamilyMatrix = buildFamilyMatrix(nFamily,eFamily,Parent);
