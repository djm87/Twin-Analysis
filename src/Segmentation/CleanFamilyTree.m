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
    G_clust.Nodes.SFAV=zeros(G_clust.numnodes,1);
    G_clust.Nodes.SFAVR=zeros(G_clust.numnodes,1);
    G_clust.Nodes.EffSF=zeros(G_clust.numnodes,1);
    exflagGroup=zeros(length(groups),1);
    progress(0,length(groups));
%     groups=588
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
%         nClustSchmidVariant = zeros(numnodes(G_clust_sub),1,'int8'); 
%         nClustK1NAng = zeros(numnodes(G_clust_sub),1);
        nSFAV = zeros(numnodes(G_clust_sub),1,'int8'); 
        nSFAVR = zeros(numnodes(G_clust_sub),1,'int8'); 
        nEffSF = zeros(numnodes(G_clust_sub),1); 
        
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
                G_clust.Nodes.Generation(nClustgroupId)=-1;
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
                progress(i,length(groups));
%                 fprintf('group:%d, progress: %d/%d\n',group,i,length(groups))
                ismg = ismultigraph(G_Family_sub);
                if ~ismg & opt.familyTree.useMinSpanTree
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
                    nParentFamilyId=zeros(max(nFamily),1);
                    [nGeneration,nType,nParentFamilyId,exflag] = fillGenerations(...
                        FamilyMatrixType,nGeneration,nType,nParentFamilyId,[],0,exflag,opt);

                    exflagGroup(i)=exflag;

                    %Store the generation and type data in G_clust for further
                    %processing and plotting
                    if exflag==0
                        
%                         [nSchmidVariant,nK1NAng]=GetSchmidVariants(G_Family_sub,nParentFamilyId,nType,opt);
%                         [G_Family,edgeList] = GetSchmidRelative(G_Family_sub,1,G_Family_sub.Nodes.meanOri,eType,opt)
                        for j=1:max(nFamily)
                            lInd=nClustFamily==j;
                            nClustType(lInd)=nType(j);
                            nClustParentFamilyId(lInd)=nParentFamilyId(j);
                            nClustGeneration(lInd)=nGeneration(j);
                            tmppair=[nParentFamilyId(j),j];
                            tmpEId=all(G_Family_sub.Edges.pairs==tmppair,2) | all(G_Family_sub.Edges.pairs==fliplr(tmppair),2);
                            isParent=G_Family_sub.Edges.pairs(tmpEId,:)==nParentFamilyId(j);
                            if ~isempty(isParent) && any(isParent)
                                nSFAV(lInd)=G_Family_sub.Edges.SFAV(tmpEId,isParent);
                                nSFAVR(lInd)=G_Family_sub.Edges.SFAVR(tmpEId,~isParent);
                                nEffSF(lInd)=G_Family_sub.Edges.sigma13(tmpEId,isParent);
                            end

                            %                             nClustK1NAng(lInd)=nK1NAng(j);
                        end

                        G_clust.Nodes.type(nClustgroupId)=nClustType;
                        G_clust.Nodes.Generation(nClustgroupId)=nClustGeneration;
                        G_clust.Nodes.ParentFamilyId(nClustgroupId)=nClustParentFamilyId;
                        G_clust.Nodes.nSFAV(nClustgroupId)=nSFAV;
                        G_clust.Nodes.nSFAVR(nClustgroupId)=nSFAVR;
                        G_clust.Nodes.nEffSF(nClustgroupId)=nEffSF;
                        %                         G_clust.Nodes.K1NAng(nClustgroupId)=nClustK1NAng;
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
%     G.Nodes.SchmidVariant=G_clust.Nodes.SchmidVariant;
    G.Nodes.nSFAV=G_clust.Nodes.nSFAV;
    G.Nodes.nSFAVR=G_clust.Nodes.nSFAVR;
    G.Nodes.nEffSF=G_clust.Nodes.nEffSF;
%     G.Nodes.K1NAng=G_clust.Nodes.K1NAng;
    G.Nodes.computeFamily = G_clust.Nodes.computeFamily;
    progress(length(groups),length(groups));

end
function [G_Family_sub] = reduceGraph(G_Family_sub,opt)
    %Create a tmp family graph for constructing the base valid Family graph
    if opt.familyTree.removeNSBE
        %Remove edges that don't share some boundary
        eMask=ones(G_Family_sub.numedges,1,'logical');
        eMask(G_Family_sub.Edges.FSgB~=0)=false;
        

        G_Family_sub_tmp=rmedge(G_Family_sub,find(eMask)); 
        nIsDetatched=indegree(G_Family_sub_tmp)==0 &...
                outdegree(G_Family_sub_tmp)==0;
        if any(nIsDetatched)   
            if  opt.familyTree.addNSBE_NodeRemoved
                %Get the nId of the unattached node
                nId=find(nIsDetatched);
            elseif opt.familyTree.addNSBE_FavorableParent
                %Add back in nodes if their orientations are among the largest volume fraction in the 
                %initial microstructure compared to other nodes in the grain
                nId=find(nIsDetatched & G_Family_sub_tmp.Nodes.FVol/max(G_Family_sub_tmp.Nodes.FVol)>0.1);
            end
            %Keep the edges associated with the removed node
            pairs=G_Family_sub_tmp.Edges.EndNodes;
            for i=1:length(nId)
                eMask(any(nId(i)==pairs,2))=false;
            end
            
%             [eid] = inedges(G_Family_sub,nId);    
%             eMask(eid)=false;
%             [eid] = outedges(G_Family_sub,nId);   
%             eMask(eid)=false;
            G_Family_sub_tmp=rmedge(G_Family_sub,find(eMask)); 
        end
    else
        G_Family_sub_tmp=G_Family_sub;
    end
    
    %Update pairs with sub graph node IDs 
    G_Family_sub_tmp.Edges.pairs=G_Family_sub_tmp.Edges.EndNodes;
    
%     figure;p=plot(G_Family_sub_tmp);


    root = findroot(G_Family_sub_tmp,opt);

    
    %Get the shortest path directed tree
%     gen=zeros(G_Family_sub_tmp.numnodes,1);
%     eLen = edgeLength(G_Family_sub_tmp,gen,opt);
%     G_Family_sub_tmp.Edges.Weight=eLen;
%     G_Family_sub_tmp.Edges.Id=[1:G_Family_sub_tmp.numedges]';
%     [TR,~,E] = shortestpathtree(G_Family_sub_tmp,root);%,'Method','unweighted');
%     
%     %Find the nodes that are no longer attached to the graph
%     nNotRemoved=unique(TR.Edges.EndNodes);
%     eMask=ones(G_Family_sub_tmp.numedges,1,'logical');
%     pairs=G_Family_sub_tmp.Edges.EndNodes;
%     for i=1:G_Family_sub_tmp.numnodes
%         if ~any(nNotRemoved==i)
%             eMask(any(i==pairs,2))=false;
%         end
%     end
%     eMask(TR.Edges.Id)=false;
%     G_Family_sub_tmp=rmedge(G_Family_sub_tmp,find(eMask)); 

%     figure;p = plot(G_Family_sub_tmp);
%     highlight(p,TR,'EdgeColor','r')
    [G_Family_sub] = minSpanFamilyTree(G_Family_sub_tmp,root,opt);
%     [T,pred,score,exflag,G_undir,gen] = minspangraph(G_Family_sub_tmp,root,opt);
%     figure;p = plot(G_undir);
%     highlight(p,T,'EdgeColor','r')
    
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
%     eMask=ones(G_Family_sub_tmp.numedges,1,'logical');
%     eMask(T.Edges.Id)=false;
%     G_Family_sub=rmedge(G_Family_sub_tmp,find(eMask)); 
%     
%     %Flip edges that so the predessesor (pred) from minspangraph is represented by the directional
%     %graph
%     eFlip = edgesToFlip(G_Family_sub,T,pred,root);
%     G_Family_sub=flipedge(G_Family_sub,eFlip); 
    
%     [TR,~,E] = shortestpathtree(G_Family_sub,root);%,'Method','unweighted');
    
%     figure;p = plot(G_Family_sub);
%     highlight(p,TR,'EdgeColor','r')
    
    
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
return
end
function root = findroot(G_Family_sub,opt)
   
    if any(G_Family_sub.Nodes.root)
        root=G_Family_sub.Nodes.Family(G_Family_sub.Nodes.root);
    else
        %Find the candidate parents using the outcloseness metric
        gen=zeros(G_Family_sub.numnodes,1);
        eLen = edgeLength(G_Family_sub,gen,opt);
        outCloseness=centrality(G_Family_sub,'outcloseness','cost',eLen);

        if opt.rootFamily.weightByFVOL
            outCloseness=outCloseness.*G_Family_sub.Nodes.FVol;
        end

        isCandidate=outCloseness>0;
        numCandidates=sum(isCandidate);
        if numCandidates==1
            root=find(isCandidate);
        else
            [~,ind] = sort(outCloseness(isCandidate));
            rootCandidates=find(isCandidate);
            rootCandidates=rootCandidates(ind);
            root=rootCandidates(end);
        end
    end
end
function eLen = edgeLength(G,gen,opt)
%EdgeLength computes the distance to traverse a twin edge relationship.
    %pull out the voting quantities
    pairs=G.Edges.pairs;
    FSgB=G.Edges.FSgB;
    FSgB(FSgB==0)=1;
    FgBL=G.Nodes.FgBL;
    EFFSF=G.Edges.EffSF(:,1);
%     EFFSF(EFFSF>0)=0;%Parent is positive so if not a parent by SF set the maximum distance
    FArea=G.Nodes.FArea;
    FVOL=G.Nodes.FVol;
    eType=G.Edges.meanType;
    if G.numedges==1
        FPArea=sum(FArea(pairs));
        FPVOL=sum(FVOL(pairs));
    else
        FPArea=sum(FArea(pairs),2);
        FPVOL=sum(FVOL(pairs),2);
    end
    
    %Create the votes
    eLen=ones(G.numedges,1);
    for i=1:size(pairs,1)
        w=opt.twin{eType(i)}.voteWeights;
        genp=gen(pairs(i,2));
        
        %Only consider texture weighting if the generation is root
        if genp==0
            eLen(i)=(1/sum(w(1:4))*(w(1)*FSgB(i)./FgBL(pairs(i,2)) +...
                 w(2)*(EFFSF(i)+0.5)/1.0 +... 
                 w(3)*FArea(pairs(i,1))./FPArea(i) +...
                 w(4)*FVOL(pairs(i,1))./FPVOL(i)))^-1; 
        else
            eLen(i)=(1 + opt.familyTree.genWeight*(genp - 1)).*...
                (1/sum(w(1:3))*(w(1)*FSgB(i)./FgBL(pairs(i,2)) +...
                 w(2)*(EFFSF(i)+0.5) +... 
                 w(3)*FArea(pairs(i,1))./FPArea(i)))^-1; 
        end
    end
    
end
function eLen = edgeLength2(G,eLen,gen,ploc,cloc,pairs,nType,opt)
%EdgeLength computes the distance to traverse a twin edge relationship.
    %Get whether the predecessor node is first or second column

    %pull out the voting quantities
    FSgB=G.Edges.FSgB;
    FSgB(FSgB==0)=1;
    FgBL=G.Nodes.FgBL;
    EFFSF=G.Edges.EffSF;
%     EFFSF(EFFSF>0)=0;%Parent is positive so if not a parent by SF set the maximum distance
    FArea=G.Nodes.FArea;
    FVOL=G.Nodes.FVol;
    eType=G.Edges.meanType;
    if G.numedges==1
        FPArea=sum(FArea(pairs));
        FPVOL=sum(FVOL(pairs));
    else
        FPArea=sum(FArea(pairs),2);
        FPVOL=sum(FVOL(pairs),2);
    end
    
    %Create the votes
    ind=find(cloc~=0 & eLen==inf);
    for i=ind'
        w=opt.twin{eType(i)}.voteWeights;
        genp=gen(pairs(i,ploc(i)));
        nTypep=nType(pairs(i,ploc(i)));
        
        if nTypep==eType(i) && opt.familyTree.noSameType
            eLen(i)=100;
        else
            %Only consider texture weighting if the generation is root
            if genp==0
                eLen(i)=(1/sum(w(1:4))*(w(1)*FSgB(i)./FgBL(pairs(i,cloc(i))) +...
                     w(2)*(EFFSF(i,ploc(i))+0.5)/1.0 +... 
                     w(3)*FArea(pairs(i,ploc(i)))./FPArea(i) +...
                     w(4)*FVOL(pairs(i,ploc(i)))./FPVOL(i)))^-1; 
            else
                eLen(i)=(1 + opt.familyTree.genWeight*(genp)).*...
                    (1/sum(w(1:3))*(w(1)*FSgB(i)./FgBL(pairs(i,cloc(i))) +...
                     w(2)*(EFFSF(i,ploc(i))+0.5) +... 
                     w(3)*FArea(pairs(i,ploc(i)))./FPArea(i)))^-1; 
            end
        end
    end
    
end
function [G_Family_sub] = minSpanFamilyTree(G_Family_sub,root,opt)
%This is modified version of prim algorithm for determining the minimum spanning tree
%The main modifications center around introducing logical related to twin formation

    %Make an undirected graph based on the nodes and edges of G_Family_sub
    G_Family_sub.Edges.Id=[1:G_Family_sub.numedges]';
    
%     FamilyMat=full(adjacency(G_Family_sub))

    
    G_undir=graph(G_Family_sub.Edges,G_Family_sub.Nodes);
    eType=G_undir.Edges.meanType;
    
    
%         figure
%     p=plot(G_Family_sub);
% %     layout(p,'force','Iterations',30)
%     p.EdgeFontSize=13;p.NodeFontSize=14;p.ArrowSize=15;p.EdgeColor=[0,0,0];
%     p.NodeColor=[1,0,0];p.MarkerSize=8;
%     pair1=G_Family_sub.Edges.EndNodes(:,1);
%     pair2=G_Family_sub.Edges.EndNodes(:,2);
%     npairs=length(pair1);
%     labeledge(p,pair1(G_undir.Edges.Id),pair2(G_undir.Edges.Id),1:npairs); 
%     
%     
    %Compute the laplacian (number of edge relationships for a node)
    L = laplacian(G_undir);
    Ldiag=diag(L);
    
    %Compute the adjacency matrix for the for indexing edges based on nodes
    %And store edge id in the matrix
    A = adjacency(G_undir);
    pairs=G_undir.Edges.EndNodes;
    for i=1:size(pairs,1)
        A(pairs(i,1),pairs(i,2))=i;
        A(pairs(i,2),pairs(i,1))=i;
    end
    
    %Initialize arrays 
    Q=zeros(G_undir.numnodes,1,'logical'); %nodes that are a part of the MST
    eLen=inf(G_undir.numedges,1); %Weights matrix
    eId=1:G_undir.numedges; %Edge id
    nType=zeros(G_undir.numnodes,1);
    gen=-ones(G_undir.numnodes,1); %Generation of node 
    pred=-ones(G_undir.numnodes,1); %Predecessor of node
    eRemove=ones(G_undir.numedges,1,'logical');
    
    ploc=zeros(G_undir.numedges,1);
    cloc=zeros(G_undir.numedges,1);
    %Handle nodes that are not connected to the graph
    Q(Ldiag==0)=true; %nodes with no edge to them
    
    %Set the root of the MST
    newNode=root;
    gen(root)=0;
    pred(root)=0;
    Q(root)=true;
    
    %Try merging the root with similar families 
%     angle(G_undir.Nodes.meanOri,G_undir.Nodes.meanOri(9))/degree
    
    %Start loop to build graph
    nNodes=length(Q)-sum(Q);
    for i=1:nNodes
        %Get edges connected to the new node
        eNew=A(newNode,:);
        
        %Build pairs to try to add
        eNew=nonzeros(eNew);
        pairsAdd=pairs(eNew,:);
        
        %If edge weight is not calculated then compute weights
        %any edge weight that is already calculated doesn't need
        %to be recalculated since it would result in a circular
        %relationship
        
        ploc(newNode==pairs(:,1))=1;
        cloc(newNode==pairs(:,1))=2;
        ploc(newNode==pairs(:,2))=2;
        cloc(newNode==pairs(:,2))=1;
        
        eLen = edgeLength2(G_undir,eLen,gen,ploc,cloc,pairs,nType,opt);
        
        %Only consider edges whose nodes aren't both in the MST 
        %and whose weights are computed. This takes care of 
        %circular relationships 
        E=zeros(G_undir.numedges,1,'logical');
        clocInd=find(cloc>0);
        for j=1:length(clocInd)
            E(clocInd(j)) = ~Q(pairs(clocInd(j),cloc(clocInd(j))));
        end
        E = E & eLen~=inf & eRemove;
        [~,ind]=min(eLen(E));
        eIdloc=eId(E);
        eIdadd=eIdloc(ind);

        if isempty(eIdadd)
            break;
        end
        
        %Update the nodes in the MST along with generation and 
        %predecessor info
        if gen(pairs(eIdadd,1))>-1
            newNode=pairs(eIdadd,2);
            oldNode=pairs(eIdadd,1);
        else
            newNode=pairs(eIdadd,1);
            
            oldNode=pairs(eIdadd,2);
        end
               
        %Update arrays
        Q(newNode)=true;
        gen(newNode)=gen(oldNode)+1;
        pred(newNode)=oldNode;
        nType(newNode)=eType(eIdadd);
        eRemove(eIdadd)=false;
    end
    
%     figure;p = plot(G_Family_sub);
    pair1=G_Family_sub.Edges.EndNodes(:,1);
    pair2=G_Family_sub.Edges.EndNodes(:,2);
    npairs=length(pair1);
%     labeledge(p,pair1,pair2,G_Family_sub.Edges.meanType); 
    
    %Create the MS
    G_MST=rmedge(G_undir,find(eRemove)); 
    G_Family_sub=rmedge(G_Family_sub,G_undir.Edges.Id(eRemove)); 

%         figure;p = plot(G_Family_sub);
    pair1=G_Family_sub.Edges.EndNodes(:,1);
    pair2=G_Family_sub.Edges.EndNodes(:,2);
    npairs=length(pair1);
%     labeledge(p,pair1,pair2,G_Family_sub.Edges.meanType); 
    
%     %Flip edges that so the predessesor (pred) from minspangraph is represented by the directional
%     %graph
    eFlip = edgesToFlip(G_Family_sub,G_MST,pred,root);

    G_Family_sub=flipedge(G_Family_sub,eFlip); 
%     figure;plot(G_Family_sub)
return
end
function [T,pred,score,exflag,G_undir,gen] = minspangraph(G_Family_sub,root,opt)
%minspangraph iterates over the Family graph given a root node   
    
    %Initialize the gen information with the minimum spanning unweighted graph
    G_Family_sub.Edges.Weight=ones(G_Family_sub.numedges,1);
    G_undir=graph(G_Family_sub.Edges,G_Family_sub.Nodes);
    [T,pred] = minspantree(G_undir,'Type','tree','Root',findnode(G_undir,root),'Method','sparse');
    gen=genFromPred(pred,zeros(length(pred),1),root,0); 
    
    %Compute the initial weights
    Weights=edgeLength(G_Family_sub,gen,opt);
    G_Family_sub.Edges.Weight=Weights;
    G_Family_sub.Edges.Id=[1:G_Family_sub.numedges]';
    G_undir=graph(G_Family_sub.Edges,G_Family_sub.Nodes);
%     G_undir.Nodes.FgBL=G_Family_sub.Nodes.FgBL;
   
%     figure;p = plot(G_Family_sub);

        
    cnt=1;
    while cnt< opt.familyTree.maxMinSpanTreeItter
        %Start loop
        [T,pred] = minspantree(G_undir,'Type','tree','Root',findnode(G_undir,root),'Method','sparse');
       

        gen=genFromPred(pred,zeros(length(pred),1),root,0);     
        
        figure;p = plot(G_undir,'layout','force');
        highlight(p,T,'EdgeColor','r')
        
        %Compute the new weights
        pairs=T.Edges.pairs;
        for k=1:T.numedges
           p1=pred(pairs(k,1));
           p2=pred(pairs(k,2));

           isOut=p2==pairs(k,:)|p1==pairs(k,:);
           if find(isOut)==2
             pairs(k,:)=fliplr(pairs(k,:));  
           end
        end
        T.Edges.pairs=pairs;
        eId=T.Edges.Id;
        Weights(eId)=edgeLength(T,gen,opt);
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
    
    score = sum(Weights(eId));
    if cnt== opt.familyTree.maxMinSpanTreeItter
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