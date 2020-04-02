function [G_Family,G_clust,G,exflagGroup]= CleanFamilyTree(groups,G_Family,G_clust,G,grains,mGrains,opt)   
%cleanFamilyTree automatically removes circular twin relationships and
%outputs a list of groups that need user input to resolve.

%exflagGroup=-1 single frag
%exflagGroup=0 clean family
%exflagGroup=1 no twin relationships but multiple families
%exflagGroup=2 too many parents
%exflagGroup=3 no parent but twin relationships exist
%exflagGroup=4 max number of generations hit while assigning generation
%exflagGroup=5 not all fragments were related... something is wrong

%==========================================================================
   
    %loop over groups
    exflagGroup=zeros(length(groups),1);
    i=1;
    while i<=length(groups)
        group=groups(i);
        
        %load edge and node properties for Family graph
        G_Family_sub = subgraph(G_Family,find(group==G_Family.Nodes.Group));
    
        %load edge and node properties for cluster graph
        nClustgroupId = find(group==G_clust.Nodes.Group);
        G_clust_sub = subgraph(G_clust,nClustgroupId);
        nClustFamily = G_clust_sub.Nodes.FamilyID;
        eClustFamily = G_clust_sub.Edges.FamilyID;
        eClustType = G_clust_sub.Edges.type;
        nClustGeneration = -ones(numnodes(G_clust_sub),1,'int8');
        nClustType = -ones(numnodes(G_clust_sub),1,'int8');      
        
        
        contInnerWhile=true;
        while contInnerWhile
            nFamily = G_Family_sub.Nodes.Family;
            eType = G_Family_sub.Edges.meanTypeRlx;
            Parent = G_Family_sub.Edges.Parent;
            eFamily = G_Family_sub.Edges.FamilyID;
            nGeneration = -ones(numnodes(G_Family_sub),1,'int8');
            nType = -ones(numnodes(G_Family_sub),1,'int8');   

            if isempty(eType)
                %Single fragment
                G_clust.Nodes.computeFamily(nClustgroupId)=false;
                exflagGroup(i)=-1;
                i=i+1;
                contInnerWhile=false;
            elseif all(eType==0)
                %Multiple fragments but no relationships
                exflagGroup(i)=1;
                i=i+1;
                contInnerWhile=false;
            else

                %Make Family relation matrix 
                FamilyMatrix = logical(full(adjacency(G_Family_sub)));

                %Find child of none (i.e. the parent) and exclude families with
                %no relationships
                TooManyParentsExist=sum(all(~FamilyMatrix,1) & ~all(~FamilyMatrix,2)')>1;  

                %Find child with more than one parent (circular relationship)
                CircularRelationshipExists=any(sum(FamilyMatrix,1)>1);

                if TooManyParentsExist
                    exflagGroup(i)=2;

                    if opt.debugFamilyTree
                        fprintf('More than one parent is apparent.. fix manually\n')

                        [G_clust,runCleanupAgain,i,exitCleanFamily] = ClusterEditor(group,G_clust,grains,mGrains,grains.meanOrientation,i,1,1,0,1,0); 
                        if exitCleanFamily
                            fprintf('exiting CleanFamilyTree\n')
                            break;
                        end
                    else
                       contInnerWhile=false;
                       i=i+1;
                    end

                elseif CircularRelationshipExists
                   [G_Family_sub] = fixCircularFamily(G_Family_sub,...
                       FamilyMatrix,eFamily,eClustFamily,eClustType,Parent,opt)

                else %Everything looks good to process!

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
                    if exflag==0
                        for j=1:max(nFamily)
                            lInd=nClustFamily==j;
                            nClustType(lInd)=nType(j);
                            nClustGeneration(lInd)=nGeneration(j);
                        end

                        G_clust.Nodes.type(nClustgroupId)=nClustType;
                        G_clust.Nodes.Generation(nClustgroupId)=nClustGeneration;
                        G_clust.Nodes.computeFamily(nClustgroupId)=false;
                    end

                    contInnerWhile=false;
                    i=i+1;
                    
                end
            end
%         G_clust.Nodes.nClustgroupId
        end
    end
    
    %Store in the initial graph so if G_clust is remade all information is
    %kept. Note, G_clust has all the nodes of G but only a subset of edges.
    G.Nodes.type = G_clust.Nodes.type;
    G.Nodes.Generation = G_clust.Nodes.Generation;
    G.Nodes.computeFamily = G_clust.Nodes.computeFamily;

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
function [G_Family_sub] = fixCircularFamily(G_Family_sub,FamilyMatrix,eFamily,eClustFamily,eClustType,Parent,opt)
    %In the case of circular twin relations a child has more than 
    %one parent and a choice between edge relationships must be made
    
    %Initialize edges
    rEdge=[]; 

    %For each child look for multiple parents
    cFList=find(sum(FamilyMatrix,1)>1);
    
    %Loop over children with too many parents
    for cF=cFList
        %For some child find all the parents
        pF=find(FamilyMatrix(:,cF));

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
            FRgB(kk)=G_Family_sub.Edges.FRgB(eId(kk),Parent(eId(kk),:));
            G_Family_sub.Edges.Vote(eId(kk),Parent(eId(kk),:));
            EffSF(kk)=G_Family_sub.Edges.EffSF(eId(kk),Parent(eId(kk),:));    
        end
        
        %Prefer first grains that are touching, then vote based on schmid
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
    if ~isempty(rEdge)
        %Remove edge and retry 
        G_Family_sub=rmedge(G_Family_sub,rEdge); 
    end
end