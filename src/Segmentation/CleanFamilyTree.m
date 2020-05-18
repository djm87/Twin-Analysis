function [G_Family,G_clust,G,exflagGroup]= CleanFamilyTree(groups,G_Family,G_clust,G,grains,mGrains,opt)   
%cleanFamilyTree automatically removes circular twin relationships and outputs a group list for
%cases that require user input

%exflagGroup=-1 single frag
%exflagGroup=0 clean family
%exflagGroup=1 no twin relationships but multiple families
%exflagGroup=2 too many parents
%exflagGroup=3 no parent but twin relationships exist
%exflagGroup=4 max number of generations hit while assigning generation
%exflagGroup=5 not all fragments were related... something is wrong
%exflagGroup=6 entered fix circular relationship fnc too many times
%exflagGroup=7 circular relationship isn't obvious enough... manually solve
%==========================================================================
   
    %loop over groups
    exflagGroup=zeros(length(groups),1);
%     groups=268
    for i=1:length(groups)
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
        nClustParentFamilyId = zeros(numnodes(G_clust_sub),1,'int8'); 
        nClustSchmidVariant = zeros(numnodes(G_clust_sub),1,'int8'); 
        nClustK1NAng = zeros(numnodes(G_clust_sub),1);
        
        loopcnt=1;maxloopcnt=3;
        while loopcnt<maxloopcnt
            nFamily = G_Family_sub.Nodes.Family;
            eType = G_Family_sub.Edges.meanType;
            Parent = G_Family_sub.Edges.Parent;
            eFamily = G_Family_sub.Edges.FamilyID;
            nGeneration = -ones(numnodes(G_Family_sub),1,'int8');
            nType = -ones(numnodes(G_Family_sub),1,'int8');   

            if isempty(eType)
                %Single fragment
                G_clust.Nodes.type(nClustgroupId)=0;
                G_clust.Nodes.Generation(nClustgroupId)=0;
                G_clust.Nodes.computeFamily(nClustgroupId)=false;
                G_clust.Nodes.ParentFamilyId(nClustgroupId)=0;
                G_clust.Nodes.SchmidVariant(nClustgroupId)=0;
                G_clust.Nodes.K1NAng(nClustgroupId)=0;
                exflagGroup(i)=-1;
                break
            elseif all(eType==0)
                %Multiple fragments but no relationships
                exflagGroup(i)=1;
                G_clust.Nodes.type(nClustgroupId)=0;
                G_clust.Nodes.Generation(nClustgroupId)=0;
                G_clust.Nodes.computeFamily(nClustgroupId)=false;
                G_clust.Nodes.ParentFamilyId(nClustgroupId)=0;
                G_clust.Nodes.SchmidVariant(nClustgroupId)=0;
                G_clust.Nodes.K1NAng(nClustgroupId)=0;
                break
            else
                fprintf('group:%d, progress: %d/%d\n',group,i,length(groups))
                ismg = ismultigraph(G_Family_sub);
                if ~ismg
                    [G_Family_sub] = reduceGraph(G_Family_sub,opt);
                end
                %Make Family relation matrix 
                FamilyMatrix = logical(full(adjacency(G_Family_sub)));

                %Find child of none (i.e. the parent) and exclude families with
                %no relationships
                TooManyParentsExist=sum(all(~FamilyMatrix,1) & ~all(~FamilyMatrix,2)')>1;  

                %Find child with more than one parent (circular relationship)
                CircularRelationshipExists=any(sum(FamilyMatrix,1)>1);
                
                
                if ismg
                   exflagGroup(i)=8;
                   break
                elseif TooManyParentsExist
                   exflagGroup(i)=2;
                   break
                elseif CircularRelationshipExists
                   [G_Family_sub,flagMSolv] = fixCircularFamily(G_Family_sub,...
                       FamilyMatrix,eFamily,eClustFamily,eClustType,Parent,opt);
                   
                   if flagMSolv
                        exflagGroup(i)=7;
                        break;
                   end
                   
                   loopcnt=loopcnt+1;
                   if loopcnt==maxloopcnt
                       exflagGroup(i)=6;
                   end
                else %Everything looks good to process!

                    %Add Type to the family matrix
                    FamilyMatrixType=int8(FamilyMatrix);
                    [r,c]=find(FamilyMatrix);

                    for j=1:length(r)
                        FamilyMatrixType(r(j),c(j))=eType(all(eFamily==[r(j),c(j)],2) | all(fliplr(eFamily)==[r(j),c(j)],2));
                    end

                    %Call labeling function
                    exflag=0;
                    nGeneration=-ones(max(nFamily),1,'int8');
                    nType=zeros(max(nFamily),1,'int8');
                    nParentFamilyId=zeros(max(nFamily),1,'int8');
                    [nGeneration,nType,nParentFamilyId,exflag] = fillGenerations(...
                        FamilyMatrixType,nGeneration,nType,nParentFamilyId,[],0,exflag,opt);

                    exflagGroup(i)=exflag;

                    %Store the generation and type data in G_clust for further
                    %processing and plotting
                    if exflag==0
                        
                        [nSchmidVariant,nK1NAng]=GetSchmidVariants(G_Family_sub,nParentFamilyId,nType,opt);
                        
                        for j=1:max(nFamily)
                            lInd=nClustFamily==j;
                            nClustType(lInd)=nType(j);
                            nClustParentFamilyId(lInd)=nParentFamilyId(j);
                            nClustGeneration(lInd)=nGeneration(j);
                            nClustSchmidVariant(lInd)=nSchmidVariant(j);
                            nClustK1NAng(lInd)=nK1NAng(j);
                        end

                        G_clust.Nodes.type(nClustgroupId)=nClustType;
                        G_clust.Nodes.Generation(nClustgroupId)=nClustGeneration;
                        G_clust.Nodes.ParentFamilyId(nClustgroupId)=nClustParentFamilyId;
                        G_clust.Nodes.SchmidVariant(nClustgroupId)=nClustSchmidVariant;
                        G_clust.Nodes.K1NAng(nClustgroupId)=nClustK1NAng;
                        G_clust.Nodes.computeFamily(nClustgroupId)=false;
                    end

                    break;                    
                end
            end
        end
    end
    
    %Store in the initial graph so if G_clust is remade all information is
    %kept. Note, G_clust has all the nodes of G but only a subset of edges.
    G.Nodes.type = G_clust.Nodes.type;
    G.Nodes.Generation = G_clust.Nodes.Generation;
    G.Nodes.ParentFamilyId =G_clust.Nodes.ParentFamilyId;
    G.Nodes.SchmidVariant=G_clust.Nodes.SchmidVariant;
    G.Nodes.K1NAng=G_clust.Nodes.K1NAng;
    G.Nodes.computeFamily = G_clust.Nodes.computeFamily;

end
function [G_Family_sub] = reduceGraph(G_Family_sub,opt)
    %Remove edges don't share some boundary
    eMask=ones(G_Family_sub.numedges,1,'logical');
    eMask(G_Family_sub.Edges.FSgB~=0)=false;
    G_Family_sub=rmedge(G_Family_sub,find(eMask)); 
    
    %Find the candidate parents 
    w=[1,1,0,0];
    gen=zeros(G_Family_sub.numnodes,1);
    eLen = edgeLength(G_Family_sub,G_Family_sub.Edges.EndNodes,gen,w,opt);
    outCloseness=centrality(G_Family_sub,'outcloseness','cost',eLen);
    outCloseness=outCloseness.*G_Family_sub.Nodes.FVol;
    isCandidate=outCloseness>0;
    numCandidates=sum(isCandidate);
    if numCandidates==1
        rootCandidates=find(isCandidate);
    else
        if numCandidates>opt.maxCanidateParents
            [~,ind] = sort(outCloseness);
            rootCandidates=ind(end+1-opt.maxCanidateParents:end);
            numCandidates=opt.maxCanidateParents;
        else
            [~,ind] = sort(outCloseness(isCandidate));
            rootCandidates=find(isCandidate);
            rootCandidates=rootCandidates(ind);
        end
    end
    outCloseness(rootCandidates);
    root=rootCandidates(end);
    %Get the shortest path directed tree
    w=[1,0,0,0];
    gen=zeros(G_Family_sub.numnodes,1);
    eLen = edgeLength(G_Family_sub,G_Family_sub.Edges.EndNodes,gen,w,opt);
    G_Family_sub.Edges.Weight=eLen;
    G_Family_sub.Edges.Id=[1:G_Family_sub.numedges]';
    %     G_Family_sub
    [TR,~,E] = shortestpathtree(G_Family_sub,root);%,'Method','unweighted');
    
    %Find the nodes that are no longer attached to the graph
    nNotRemoved=unique(TR.Edges.EndNodes);
    eMask=ones(G_Family_sub.numedges,1,'logical');
    pairs=G_Family_sub.Edges.EndNodes;
    for i=1:G_Family_sub.numnodes
        if ~any(nNotRemoved==i)
            eMask(any(i==pairs,2))=false;
        end
    end
    eMask(TR.Edges.Id)=false;
    G_Family_sub2=rmedge(G_Family_sub,find(eMask)); 

    figure;p = plot(G_Family_sub2);
    highlight(p,TR,'EdgeColor','r')
    
    [T,pred,score,exflag,G_undir,gen] = minspangraph(G_Family_sub2,root,opt);

    
    %Find the best performing parent
%     T={};score={};exflag={};pred={};G_undir={};gen={};
%     for i=1:numCandidates
%         [T{i},pred{i},score{i},exflag{i},G_undir{i},gen{i}] = minspangraph(G_Family_sub,rootCandidates(i),opt);
%     end
%     score=max([gen{:}])'.*outCloseness(rootCandidates(1))./outCloseness(rootCandidates).*[score{:}]';
%     [~,ind]=min(score);
%     T=T{ind};
%     G_undir=G_undir{i};
%     pred=pred{ind};
%     score=score(ind);
%     exflag=exflag{ind};
    %Remove edges that aren't a part of the min spanning graph
    eMask=ones(G_Family_sub2.numedges,1,'logical');
    eMask(T.Edges.Id)=false;
    G_Family_sub=rmedge(G_Family_sub2,find(eMask)); 
%     
%     %Flip edges that so the predessesor (pred) from minspangraph is represented by the directional
%     %graph
    eFlip = edgesToFlip(G_Family_sub,T,pred,root);
    G_Family_sub=flipedge(G_Family_sub,eFlip); 
%     G_Family_sub=TR;
% 
%     figure;
%     p = plot(G_undir);
%     highlight(p,T);
    
%     figure;p = plot(TR);
%     figure;p=plot(G_Family_sub);
%     highlight(p,G_Family_sub,'EdgeColor','r')
% 
%     figure;
%     p = plot(G_Family_sub);
end
function eLen = edgeLength(G,pairs,gen,w,opt)
    FSgB=G.Edges.FSgB;
    FSgB(FSgB==0)=1;
    FgBL=G.Nodes.FgBL;
    EFFSF=G.Edges.EffSF(:,1);
    EFFSF(EFFSF>0)=0;%Parent is positive so if not a parent by SF set the maximum distance
    FArea=G.Nodes.FArea;
    if G.numedges==1
        FPArea=sum(FArea(pairs));
    else
        FPArea=sum(FArea(pairs),2);
    end
    
    FVOL=G.Nodes.FVol;
    if max(gen)>opt.maxGen
        %Set a bad score
        eLen=(w(1)*FSgB./FgBL(pairs(:,2)) + w(2)*abs(EFFSF)/0.5 + w(3)*FArea(pairs(:,1))./FPArea + w(4)*FVOL(pairs(:,1))./FVOL(pairs(:,2))).^-1;
        eLen=2*eLen;
    elseif all(gen(pairs(:,2)))==0
        eLen=(w(1)*FSgB./FgBL(pairs(:,2)) + w(2)*abs(EFFSF)/0.5 + w(3)*FArea(pairs(:,1))./FPArea + w(4)*FVOL(pairs(:,1))./FVOL(pairs(:,2))).^-1; 
    else
        eLen=0.8*gen(pairs(:,2)).*(w(1)*FSgB./FgBL(pairs(:,2)) + w(2)*abs(EFFSF)/0.5 + w(3)*FArea(pairs(:,1))./FPArea + w(4)*FVOL(pairs(:,1))./FVOL(pairs(:,2))).^-1; 
    end
    end
function [T,pred,score,exflag,G_undir,gen] = minspangraph(G_Family_sub,root,opt)
%minspangraph iterates over the Family graph given a root node   
    w=[1,1,1,0];
    gen=zeros(G_Family_sub.numnodes,1);
    Weights=edgeLength(G_Family_sub,G_Family_sub.Edges.EndNodes,gen,w,opt);
    G_Family_sub.Edges.Weight=Weights;
    G_Family_sub.Edges.Id=[1:G_Family_sub.numedges]';
    G_undir=graph(G_Family_sub.Edges,G_Family_sub.Nodes);
%     G_undir.Nodes.FgBL=G_Family_sub.Nodes.FgBL;

    cnt=1;
    while cnt<opt.spanTreeItter
        %Start loop
        [T,pred] = minspantree(G_undir,'Type','tree','Root',findnode(G_undir,root));
       

        gen=genFromPred(pred,zeros(length(pred),1),root,0);     
        
        %Compute the new weights
        pairs=T.Edges.EndNodes;
        for k=1:T.numedges
           p1=pred(pairs(k,1));
           p2=pred(pairs(k,2));

           isOut=p2==pairs(k,:)|p1==pairs(k,:);
           if find(isOut)==2
             pairs(k,:)=fliplr(pairs(k,:));  
           end
        end
        w=[1,0,0,0];
        eId=T.Edges.Id;
        Weights(eId)=edgeLength(T,pairs,gen,w,opt);
%         [pairs,Weights(eId)]
        if cnt~=1
            if isomorphism(T,Told)
                break;
            end
        end   
        
        %Update values for next iteration
        Told=T;
        G_Family_sub.Edges.Weight=Weights;
        G_undir=graph(G_Family_sub.Edges,G_Family_sub.Nodes);
        G_undir.Nodes.FgBL=G_Family_sub.Nodes.FgBL;
        cnt=cnt+1;
    end
    
    score = sum(Weights(eId).^-1);
    if cnt== opt.spanTreeItter
       exflag=1; 
    else
        exflag=0;
    end
end
function [gen]=genFromPred(pred,gen,parent,lvl)    
    lvl=lvl+1;
    child=find(pred==parent);
    gen(child)=lvl;
    for i=1:length(child)
        gen=genFromPred(pred,gen,child(i),lvl);  
    end
end
function [eFlip] = edgesToFlip(G_Family_sub,T,pred,root)
%edgesToFlip determines the edge direction from the min span tree and determines if the edges in the 
%directional graph need to be flipped    
    toFlip=zeros(G_Family_sub.numedges,1,'logical');
    allPairs=G_Family_sub.Edges.EndNodes; 
    pairs=T.Edges.EndNodes;
    
    for k=1:T.numedges
       p1=pred(pairs(k,1));
       p2=pred(pairs(k,2));
       if p2==0 || p1==0
           isOut=root==pairs(k,:);
       else
           isOut=p2==pairs(k,:)|p1==pairs(k,:);
       end
       
       allId=find(all(pairs(k,:)==allPairs,2)| all(fliplr(pairs(k,:))==allPairs,2));
       toFlip(allId) = pairs(k,isOut)==allPairs(allId,2);
    end

    eFlip=find(toFlip);
end
function [nGeneration,nType,nParentFamilyId,exflag] = fillGenerations(FamilyMatrix,...
                                        nGeneration,nType,nParentFamilyId,genId,gen,exflag,opt)
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
            nParentFamilyId(genId)=0;
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
            nParentFamilyId(:)=0;
            return
        end
    end
    if gen>opt.maxGen
        %max number of generations hit
        exflag=4;
        
        %Reset nGeneration so we know to debug.
        nGeneration(:)=-1;
        nType(:)=0;
        nParentFamilyId(:)=0;
        return
    else
        for i=1:length(genId)
            genId1 = find(FamilyMatrix(genId(i),:)~=0);
            if ~isempty(genId1)
                nGeneration(genId1)=gen;
                nType(genId1)=FamilyMatrix(genId(i),genId1);
                nParentFamilyId(genId1)=genId(i);
                if any(nGeneration==-1)
                    [nGeneration,nType,nParentFamilyId,exflag] = fillGenerations(FamilyMatrix,...
                        nGeneration,nType,nParentFamilyId,genId1,gen+1,exflag,opt);
                end
            end
        end
        if any(nGeneration==-1 & nType~=opt.twinUnknown) && gen==1
            if exflag~=4
               exflag=5;
            end
            nGeneration(:)=-1;
            nType(:)=0;
            nParentFamilyId(:)=0;
        end
    end
end
function [G_Family_sub,flagMSolv] = fixCircularFamily(G_Family_sub,FamilyMatrix,eFamily,eClustFamily,eClustType,Parent,opt)
    %In the case of circular twin relations a child has more than 
    %one parent and a choice between edge relationships must be made
    
    %Initialize edges
    rEdge=[]; 
    
    %For each child look for multiple parents
    cFList=find(sum(FamilyMatrix,1)>1);
    
    %Loop over children with too many parents
    flagMSolv=false;
    for cF=cFList
        %For some child find all the parents
        pF=find(FamilyMatrix(:,cF));

        %Find the right edge for the parent/child relationships
        eId=zeros(length(pF),1);
        FSgB=zeros(length(pF),1);
        EffSF=zeros(length(pF),1);
        localEdges=zeros(length(pF),1);
        
        for kk=1:length(pF)
            PotentialEdges=all(eClustFamily==[pF(kk),cF] | eClustFamily==[cF,pF(kk)],2);
            PotentialEdges(eClustType==opt.twinUnknown)=false;     
            localEdges(kk)=sum(PotentialEdges);
            eId(kk)=find(sum(or(pF(kk)==eFamily,cF==eFamily),2)==2);
            FSgB(kk)=G_Family_sub.Edges.FSgB(eId(kk));
            G_Family_sub.Edges.Vote(eId(kk),Parent(eId(kk),:));
            EffSF(kk)=G_Family_sub.Edges.EffSF(eId(kk),Parent(eId(kk),:));    
        end
        
        %Prefer first grains that are touching, then vote based on schmid
        if sum(localEdges~=0)>1 & opt.CircularRelationships.UseBoundary
            %Then all families share boundary decide based on boundary
            [FSgB_sorted,FSgB_I]=sort(FSgB);
            
            pgBdiff=abs(FSgB_sorted(end)-FSgB_sorted(end-1))/FSgB_sorted(end);
            
            if pgBdiff > opt.CircularRelationships.mingBpdiff
                rEdge=[rEdge;eId(FSgB~=FSgB_sorted(end))];
            else
                [EffSF_sorted,EffSF_I]=sort(EffSF);
                
                pEFFSFdiff=abs(EffSF_sorted(end)-EffSF_sorted(end-1))/EffSF_sorted(end);
                
                if pEFFSFdiff > opt.CircularRelationships.minEFFSFpdiff
                    rEdge=[rEdge;eId(EffSF~=EffSF_sorted(end))]; 
                else
                    flagMSolv=true;
                end
            end
            
        elseif sum(localEdges~=0)==1 & opt.CircularRelationships.UseBoundary
            %Assign the grain that shares boundary
            rEdge=[rEdge;eId(localEdges==0)];
        elseif opt.CircularRelationships.UseEFFSF
            %Use the schmid factor to figure it out
            [EffSF_sorted,EffSF_I]=sort(EffSF);

            pEFFSFdiff=abs(EffSF_sorted(end)-EffSF_sorted(end-1))/EffSF_sorted(end);

            if pEFFSFdiff > opt.CircularRelationships.minEFFSFpdiff
                rEdge=[rEdge;eId(EffSF~=EffSF_sorted(end))]; 
            else
                flagMSolv=true;
            end
        else
            flagMSolv=true;
        end
    end 
    if ~isempty(rEdge)
        %Remove edge and retry 
        G_Family_sub=rmedge(G_Family_sub,rEdge); 
    end
end