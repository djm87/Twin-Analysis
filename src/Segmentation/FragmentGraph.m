function [G,mGrains,twinGb] = FragmentGraph(grains,opt)
    %This function initializes the graph objects for grain fragments.  
    %Initialization entails creating a base graph from which the cluster
    %graph is constructed from. The starting list of "twin edges" to construct
    %the cluster graph are computed here

    %Initialize the fragment graph
    %======================================================================
    %Compute grain neighbors
    [pairs] = neighbors(grains);

    %Initialize graph
    s=pairs(:,1);
    t=pairs(:,2);
    G=graph(s,t);

    % Add node labels and other grain level information 
    nGrains=length(grains);
    G.Nodes.Id=grains.id;
    G.Nodes.centroids=grains.centroid;
    G.Nodes.Group=zeros(nGrains,1,'int32');
    G.Nodes.FamilyID=zeros(nGrains,1,'int16');
    G.Nodes.EffSF=zeros(nGrains,opt.nMori,'double');
    G.Nodes.FArea=zeros(nGrains,1,'double');
    G.Nodes.sigma13=zeros(nGrains,1,'double');
    G.Nodes.Generation = -ones(length(G.Nodes.Id),1,'int8');
    G.Nodes.ParentFamilyId =zeros(nGrains,1,'int8');
    G.Nodes.SchmidVariant=zeros(nGrains,1,'int8');
    G.Nodes.K1NAng=zeros(nGrains,1,'double');
    G.Nodes.type=zeros(nGrains,1,'int8');
    G.Nodes.typeColored = zeros(nGrains,3,'int8');
    
    %These groups store whether a cluster group or Family  has been modified 
    %since it was last computed. This reduces the amount of needless
    %computation and significantly speeds up the cluster editing process.
    G.Nodes.isNewGroup=ones(nGrains,1,'logical');
    G.Nodes.computeTree=ones(nGrains,1,'logical');

    
    % Add neighbor information
    nPairs=size(pairs,1);
    G.Edges.pairs=pairs; %these node pairs correspond to grains.id
    G.Edges.GlobalID=int32([1:nPairs]');
    G.Edges.Parent=zeros(nPairs,2,'logical');
    G.Edges.SFCalcd=zeros(nPairs,1,'logical');
    G.Edges.EffSF=zeros(nPairs,2,'double');
    G.Edges.EffSFRelative=zeros(nPairs,1,'double'); %Used?
    G.Edges.FamilyID=zeros(nPairs,2,'int16');
    G.Edges.Vote=zeros(nPairs,2,'double');
    G.Edges.FRArea=zeros(nPairs,2,'double');
    G.Edges.FRgB=zeros(nPairs,2,'double');
    G.Edges.FREffSF=zeros(nPairs,2,'double');
    G.Edges.Group=zeros(nPairs,1,'int32');
    G.Edges.meanType=zeros(nPairs,1,'int8');
    G.Edges.sigma13=zeros(nPairs,2,'double');
    G.Edges.GbInd=zeros(nPairs,1,'int32');
    G.Edges.type=zeros(nPairs,1,'int8');
    
    %Initialize G_clust quantities that can be calculated once
    %======================================================================
    %Compute the initial cluster sub graph. Edges are merged based on the
    %combine logical variable.
    [G.Edges.combine,mGrains,twinGb,G.Edges.GBType] = ...
        MergeByBoundary(G,grains,opt.gclust.beta,opt);
    [G.Edges.combineRlx,mGrainsRlx,twinGbRlx,G.Edges.GBTypeRlx] = ...
        MergeByBoundary(G,grains,opt.gclust.betaRlx,opt);

    %Get the type of misorientation (e.g. twin type 1, etc..
    mori=inv(grains(G.Edges.pairs(:,1)).meanOrientation).*...
        grains(G.Edges.pairs(:,2)).meanOrientation; 
    [~,G.Edges.meanType] = TestTwinRelationship(mori,opt.gfam.phi,opt,G.Edges.type);
    
    %Update the types to include graph segmentation angle merging 
    G = MergeLAGB(G,mori,opt);

    %Handle merging of inclusions
    if  opt.gclust.mergeFrag
        [G.Edges.typeIni,G.Edges.combine] = ...
            MergeInclusions(G.Edges.pairs,grains,G.Edges.combine,...
            opt.moriUnknown,G.Edges.meanType,G.Edges.GBType,G.Edges.GBTypeRlx);
    else
        G.Edges.typeIni=G.Edges.GBType;
    end
    
    %Ensure no zero types (This shouldn't happen if logic is airtight).
    assert(all(G.Edges.typeIni(G.Edges.combine)>0),'There is an error in type definition in FragmentGraph')
    
end

function [G] = MergeLAGB(G,mori,opt)
    toMerge=mori.angle <  opt.gclust.betaFrag;
    G.Edges.GBTypeRlx(toMerge)=opt.moriLAGB;
    G.Edges.GBType(toMerge)=opt.moriLAGB;
    G.Edges.meanType(toMerge)=opt.moriLAGB;
    G.Edges.combineRlx(toMerge)=true;
    G.Edges.combine(toMerge)=true;
end

function [type,combine] = MergeInclusions(pairs,grains,combine,moriUnknown,meanType,type,GBTypeRlx)
    %Give inclusions an edge relationship regardless of whether they end up
    %being a twin. Usually they are twins, bad indexed pixels, or an
    %artifact of a 3D microstructure.
    
    isInside = checkInside(grains, grains);
    [GrainIdInclusion,GrainIdWithInclusion] = find(isInside);

%     figure;plot(grains,grains.meanOrientation);hold on
%     text(grains,int2str(grains.id));hold off
    
    %Get pairs to consider
    epairs=[GrainIdInclusion,GrainIdWithInclusion];
    epairs=epairs(epairs(:,1)~=epairs(:,2),:);

    %Build the boundary for testing if fully included in fragment
    gBId=grains(GrainIdInclusion).boundary(grains.mineral,grains.mineral).grainId;
    gBId=unique(vertcat(gBId,fliplr(gBId)),'rows');
    for i=1:size(epairs,1)
        %Only consider true internal grains
        if sum(epairs(i,1)==gBId(:,1))==1
            eInd=find(all(epairs(i,:)==pairs,2) | all(fliplr(epairs(i,:))==pairs,2));
            isTwinGB=GBTypeRlx(eInd)~=0;
            isTwinMean=meanType(eInd)~=0;
            if isTwinGB
                type(eInd)=GBTypeRlx(eInd);
                combine(eInd)=true;
            elseif isTwinMean
                type(eInd)=meanType(eInd);
                combine(eInd)=isTwinMean;
            else
                type(eInd)=moriUnknown;
                combine(eInd)=true;
            end
        end
    end
end


