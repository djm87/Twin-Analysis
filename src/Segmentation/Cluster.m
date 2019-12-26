function [G_clust,G,mGrains] = Cluster(G,grains,opt)
    %Cluster performs the grouping based on edge definitions and user input
    %files. 
    
    % Test if all the files exist and generate if they don't
    flist={'eAddList.txt','eRemoveList.txt','eRemoveCleanup.txt','nRemoveList.txt','nAddList.txt','eRelation.txt','notParent.txt','notTwin.txt'};
    generateInputs(flist); %move to cluster
    
    eAddList=load(flist{1}); %Adds edge and tries to give it a twin relationship, otherwise type twin unknown
    eRemoveList=load(flist{2}); %Remove the edge 
    eRemoveCleanup=load(flist{3}); %Edges determined by cleanFamilyTree
    nRemoveList=load(flist{4}); %Removes all edges associated with node
    nAddList=load(flist{5}); %Tries to add all edges that have twin relationship
        
    %Assign the type using the relaxed mean misorientation tolerance.
    meanTypeRlx=G.Edges.meanTypeRlx;
    G.Edges.type=G.Edges.meanTypeRlx;
    
    %As a starting point, use tight boundary misorientation merging 
    G.Edges.combineCleaned = G.Edges.combine;
    
    %Any edges that aren't a twin should be set to be removed
    %Filter edge by mistol based with the type of meanMistolRelaxed
    G.Edges.combineCleaned(G.Edges.type==0)=false;
    
    %Try adding nodes
    pairs=G.Edges.pairs;
    for i=1:length(nAddList)
        toAdd=any(nAddList(i)==pairs,2);
        G.Edges.combineCleaned(toAdd)=G.Edges.combineRlx(toAdd);
    end
    
    %Try removing nodes
    for i=1:length(nRemoveList)
        [~,epairs]=grains(nRemoveList(i)).neighbors;
        for j=1:size(epairs,1)
            eId=find(all(epairs(j,:)==G.Edges.pairs,2) | all(fliplr(epairs(j,:))==G.Edges.pairs,2));
            isTwin=G.Edges.type(eId)~=0 & G.Edges.type(eId)~=opt.twinUnknown;
            G.Edges.combineCleaned(eId(isTwin))=false;
        end
    end
    
    %Add edges
    eAddId=cell(size(eAddList,1),1);
    for i=1:size(eAddList,1)
        eAddId{i}=find(all(eAddList(i,:)==G.Edges.pairs,2) | all(fliplr(eAddList(i,:))==G.Edges.pairs,2));
    end
    eAddId=cell2mat(eAddId);
    typeAdd=meanTypeRlx(eAddId);
    typeAdd(typeAdd==0)=opt.twinUnknown; %Unknown Type
    G.Edges.type(eAddId)=typeAdd;
    G.Edges.combineCleaned(eAddId)=true;
    
    %Update unknown types in nodes
    typeKnown=unique(G.Edges.pairs([G.Edges.type~=opt.twinUnknown,G.Edges.type~=opt.twinUnknown]));
    G.Nodes.typeUnknown(:)=true;
    G.Nodes.typeUnknown(typeKnown)=false;
  
    %Remove edges
    G.Edges.combineCleaned(eRemoveList)=false; %
    G.Edges.combineCleaned(eRemoveCleanup)=false; %

    %Merge grains based on edge combine list
    [mGrains,parentId,mTwinGb,G.Edges.GbInd] = MergeByEdge(G.Edges.pairs,G.Edges.GbInd,G.Edges.combineCleaned,grains);    
    if opt.mergeInclusionCluster
        %Add grains that are internal to a cluster
        isInside = checkInside(mGrains, mGrains);
        [GrainIdInclusion,GrainIdWithInclusion] = find(isInside);
        uGrainIdInclusion=unique(GrainIdInclusion);
        nInternal=zeros(length(grains),1,'logical');
        for i=1:length(uGrainIdInclusion)
            nInternal(uGrainIdInclusion(i)==parentId)=true;
        end
        nAddList=find(nInternal);

        %Try adding nodes using any edge that has a twin relationship with
        %Relaxed tolerance else add one edge of unknown type
        pairs=G.Edges.pairs;
        for i=1:length(nAddList)
            toAdd=any(nAddList(i)==pairs,2);
            combine=G.Edges.combineRlx(toAdd);
            if ~any(combine)
               %Then use the first unknown twin relationship to merge into
               %cluster. Note this addition will go away if nAddList doesn't
               %contain the same node in futures calls to cluster or the 
               %first edge is in edge remove list.
               combine(1)=true; 
               G.Edges.type(toAdd(1))=opt.twinUnknown;
            end
            G.Edges.combineCleaned(toAdd)=combine;
        end    

        %Try removing nodes so no internal grains are added that we specify to
        %remove
        for i=1:length(nRemoveList)
            [~,epairs]=grains(nRemoveList(i)).neighbors;
            for j=1:size(epairs,1)
                eId=find(all(epairs(j,:)==G.Edges.pairs,2) | all(fliplr(epairs(j,:))==G.Edges.pairs,2));
                isTwin=G.Edges.type(eId)~=0 & G.Edges.type(eId)~=length(twin);
                G.Edges.combineCleaned(eId(isTwin))=false;
            end
        end

        %Remove edges similarly to ensure nothing is added that we said not to
        %add
        G.Edges.combineCleaned(eRemoveList)=false; %
        G.Edges.combineCleaned(eRemoveCleanup)=false;

        %Merge grains based on edge combine list
        [mGrains,parentId,mTwinGb,G.Edges.GbInd] = MergeByEdge(G.Edges.pairs,G.Edges.GbInd,G.Edges.combineCleaned,grains);    
    end
    %Remove edges to make reduced graph over the clusters
    G_clust=rmedge(G,G.Edges.pairs(~G.Edges.combineCleaned,1),...
        G.Edges.pairs(~G.Edges.combineCleaned,2));
    
    %Assign group Ids using mGrains id
    [G_clust.Nodes.Group,G_clust.Edges.Group] = GroupsFromGraph(G_clust.Edges.pairs,length(grains),length(mGrains),parentId);
    
    %Transfer still valid Family and family vote info to G_clust.. build
    %list of groups that need to be computed
    [groupList] = GetGroupsToCompute(G_clust,G);
    
    %Update G with groups in G_clust
    G.Nodes.Group(G_clust.Nodes.Id)=G_clust.Nodes.Group;
    G.Edges.Group(G.Edges.combineCleaned)=G_clust.Edges.Group;
    
    %Update unknown nodes
    nId_edges=unique(G_clust.Edges.pairs);
    G_clust.Nodes.typeUnknown(:)=true; 
    G_clust.Nodes.typeUnknown(nId_edges)=false;

    % Identify Families in grain clusters 
    [G_clust,G]= AssignFamilyIDs(G_clust,G,groupList,grains,mGrains,opt);
    
    %Compute Schmid info for twin/parents in clustered grains
    %This computes Schmid factor for twin/parent identification
    [G_clust,G]= GetSchmidRelative(G_clust,G,grains,mGrains,opt);

    %Perform family votes
    [G_clust,G] = FamilyVotes(G_clust,G,groupList,grains,mGrains,opt);    

    %Plot the results
    if opt.plot.do
        %Remove grains that haven't twinned (better plot performance)
        if opt.plot.ClusterOnly
            tokeep=zeros(length(G_clust.Nodes.Id),1,'logical');
            tokeep(unique(G_clust.Edges.pairs))=true;
            G_clust_clean=rmnode(G_clust,find(~tokeep));
        else
            G_clust_clean=G_clust;
        end
        
        %Plot edge labeled graph
        labelNodes=false;labelEdges=opt.plot.labelEdges;plotG=true;legendOn=false;
        fhandle = plotGraph(grains,mGrains,G_clust_clean,...
            grains.meanOrientation,G_clust_clean.Nodes.Id,...
            labelNodes,labelEdges,legendOn,plotG,[]);
        
        %Plot node labeled graph
        labelNodes=true;labelEdges=false;plotG=false;legendOn=false;
        fhandle = plotGraph(grains,mGrains,G_clust_clean,...
            grains.meanOrientation,G_clust_clean.Nodes.Id,...
            labelNodes,labelEdges,legendOn,plotG,[]);
        
    end
end
function [groupList] = GetGroupsToCompute(G_clust,G)
    %Get group ids - optimized 12/24/19
    groupCompute=zeros(max(G_clust.Nodes.Group),1,'logical');
    for i=1:max(G_clust.Nodes.Group)
        %Extract group to local variables for readability
        eID=G_clust.Edges.GlobalID(i==G_clust.Edges.Group);
        nID = find((i==G_clust.Nodes.Group));
        G_nID = find(G.Nodes.Group(nID(1))==G.Nodes.Group);
        G_eId = find(G.Nodes.Group(nID(1))==G.Edges.Group);
               
        if length(G_nID) ~= length(nID) || length(G_eId) ~= length(eID) 
            groupCompute(i)=true;
        elseif ~all(nID==G_nID) || ~all(eID==G_eId) %Because sorted
            groupCompute(i)=true;           
        end
    end
    groupList=find(groupCompute);
end
