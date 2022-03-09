function [G,mGrains,twinGb] = InitialGraph(ebsd,grains,opt)
    %This function initializes graph objects for grains.  
    %Initialization entails creating a base graph from which the cluster
    %graph is constructed from.

    %Compute grain neighbors
    pairs = neighbors(grains);

    %Initialize graph
    s=pairs(:,1);
    t=pairs(:,2);
    G=graph(s,t);

    % Add node labels and other grain level information 
    G.Nodes.Id=grains.id;
    G.Nodes.centroids=grains.centroid;
    G.Nodes.Group=zeros(length(grains),1,'int32');
    G.Nodes.FamilyID=zeros(length(grains),1,'int16');
    G.Nodes.EffSF=zeros(length(grains),opt.nTwin,'double');
    G.Nodes.FArea=zeros(length(grains),1,'double');
%     G.Nodes.typeUnknown=zeros(length(grains),1,'logical');
    G.Nodes.sigma13=zeros(length(grains),1,'double');
    G.Nodes.Generation = -ones(length(G.Nodes.Id),1,'int8');
    G.Nodes.ParentFamilyId =zeros(length(grains),1,'int8');
    G.Nodes.SchmidVariant=zeros(length(grains),1,'int8');
    G.Nodes.K1NAng=zeros(length(grains),1,'double');
    G.Nodes.type=zeros(length(grains),1,'int8');
    G.Nodes.typeColored = zeros(length(G.Nodes.Id),3,'int8');
    
    %These groups say where a cluster group or Family  has been modified 
    %since it was last computed. This reduces the amount of needless
    %computation and significantly speeds up the cluster editing process.
    G.Nodes.isNewGroup=ones(length(G.Nodes.Id),1,'logical');
    G.Nodes.computeFamily=ones(length(G.Nodes.Id),1,'logical');
    
    % Add intergranular information
    G.Edges.pairs=pairs; %these node pairs correspond to grains.id
    G.Edges.GlobalID=int32([1:size(pairs,1)]');
    G.Edges.Parent=zeros(size(pairs,1),2,'logical');
    G.Edges.SFCalcd=zeros(size(pairs,1),1,'logical');
    G.Edges.EffSF=zeros(size(pairs,1),2,'double');
    G.Edges.EffSFRelative=zeros(size(pairs,1),1,'double'); %Used?
    G.Edges.FamilyID=zeros(size(pairs,1),2,'int16');
    G.Edges.Vote=zeros(size(pairs,1),2,'double');
    G.Edges.FRArea=zeros(size(pairs,1),2,'double');
    G.Edges.FRgB=zeros(size(pairs,1),2,'double');
    G.Edges.FREffSF=zeros(size(pairs,1),2,'double');
    G.Edges.Group=zeros(size(pairs,1),1,'int32');
    G.Edges.meanType=zeros(size(pairs,1),1,'int8');
    G.Edges.sigma13=zeros(size(pairs,1),2,'double');
    G.Edges.GbInd=zeros(size(pairs,1),1,'int32');
    G.Edges.type=zeros(size(pairs,1),1,'int8');
    
    %Compute the initial cluster sub graph. Edges are merged based on the
    %combine logical variable.
    misTol=zeros(opt.nTwin,1);
    misTolRlx=zeros(opt.nTwin,1);
    for i=1:opt.nTwin
        misTol(i)=opt.twin{i}.tol.misGb;
        misTolRlx(i)=opt.twin{i}.tol.misGbRlx;
    end
    [G.Edges.combine,mGrains,twinGb,G.Edges.GBType] = ...
        MergeByBoundary(G,grains,misTol,opt);
    [G.Edges.combineRlx,mGrainsRlx,twinGbRlx,G.Edges.GBTypeRlx] = ...
        MergeByBoundary(G,grains,misTolRlx,opt);

    %Get the type of misorientation (e.g. twin type 1, etc..
    mori=inv(grains(G.Edges.pairs(:,1)).meanOrientation).*...
        grains(G.Edges.pairs(:,2)).meanOrientation; 
    meanMisTol=zeros(opt.nTwin,1);
    for i=1:opt.nTwin
        meanMisTol(i)=opt.twin{i}.tol.misMean;
    end
    [~,G.Edges.meanType] = TestTwinRelationship(mori,meanMisTol,opt,G.Edges.type);
    
    %Update the types to include graph segmentation angle merging 
    toMerge=mori.angle <  opt.grain_recon.segAngle;
    G.Edges.GBTypeRlx(toMerge)=opt.grain_recon.segType;
    G.Edges.GBType(toMerge)=opt.grain_recon.segType;
    G.Edges.meanType(toMerge)=opt.grain_recon.segType;
    G.Edges.combineRlx(toMerge)=true;
    G.Edges.combine(toMerge)=true;

    %Handle merging of inclusions
    if  opt.mergeByGeometry.mergeFrag
        [G.Edges.typeIni,G.Edges.combine] = ...
            AutoMergeInclusions(G.Edges.pairs,grains,G.Edges.combine,...
            opt.twinUnknown,G.Edges.meanType,G.Edges.GBType,G.Edges.GBTypeRlx);
    else
        G.Edges.typeIni=G.Edges.GBType;
    end
    
    %Ensure no zero types (DS: This shouldn't happen if logic is airtight).
    assert(all(G.Edges.typeIni(G.Edges.combine)>0),'There is an error in type definition in initialGraph')
    
    end

function [type,combine] = AutoMergeInclusions(pairs,grains,combine,unknownType,meanType,type,GBTypeRlx)
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
                type(eInd)=unknownType;
                combine(eInd)=true;
            end
        end
    end
end


