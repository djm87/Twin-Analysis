function [G,mGrains,mGrainsRlx] = InitialGraph(ebsd,grains,opt)
    %This function initializes graph objects for grains.  
    %Initialization entails computation of all local grain properties such as
    %orientation, centroid, aspectRatio, etc.. and is stored in edge pair
    %(intergranular) or nodes (grains). 

    %Compute grain neighbors
    [counts,pairs] = neighbors(grains);

    %Initialize graph
    s=pairs(:,1);
    t=pairs(:,2);
    G=graph(s,t);

    % Add node labels and other grain level information (many of these
    % should be removed since the goal after segmentation is just to use
    % the grains produced by MTEX.
    G.Nodes.Id=grains.id;
    G.Nodes.centroids=grains.centroid;
    G.Nodes.Group=zeros(length(grains),1,'int32');
    G.Nodes.FamilyID=zeros(length(grains),1,'int16');
    G.Nodes.EffSF=zeros(length(grains),opt.nTwin,'double');
    G.Nodes.FArea=zeros(length(grains),1,'double');
    G.Nodes.typeUnknown=zeros(length(grains),1,'logical');
    G.Nodes.sigma13=zeros(length(grains),1,'double');

    % Add intergranular information
    G.Edges.pairs=pairs; %this contains the grain ids corresponding to grains.id
    G.Edges.GlobalID=int32([1:size(pairs,1)]');
    G.Edges.Parent=zeros(size(pairs,1),2,'logical');
    G.Edges.SFCalcd=zeros(size(pairs,1),1,'logical');
    G.Edges.EffSF=zeros(size(pairs,1),2,'double');
    G.Edges.EffSFRelative=zeros(size(pairs,1),1,'double');
    G.Edges.FamilyID=zeros(size(pairs,1),2,'int16');
    G.Edges.Vote=zeros(size(pairs,1),2,'double');
    G.Edges.FRArea=zeros(size(pairs,1),2,'double');
    G.Edges.FRgB=zeros(size(pairs,1),2,'double');
    G.Edges.FREffSF=zeros(size(pairs,1),2,'double');
    G.Edges.Group=zeros(size(pairs,1),1,'int32');
    G.Edges.meanType=zeros(size(pairs,1),1,'int8');
    G.Edges.meanTypeRlx=zeros(size(pairs,1),1,'int8');
    G.Edges.sigma13=zeros(size(pairs,1),2,'double');
    G.Edges.GbInd=zeros(size(pairs,1),1,'int32');
    G.Edges.type=zeros(size(pairs,1),1,'int8');
    G.Nodes.Generation = -ones(length(G.Nodes.Id),1,'int8');
    G.Nodes.typeColored = zeros(length(G.Nodes.Id),3,'int8');
    G.Nodes.isNewGroup=ones(length(G.Nodes.Id),1,'logical');
    G.Nodes.computeFamily=ones(length(G.Nodes.Id),1,'logical');
    
    tol=zeros(opt.nTwin,1);
    tolRlx=zeros(opt.nTwin,1);
    for i=1:opt.nTwin
        tol(i)=opt.twin{i}.tol.misGb;
        tolRlx(i)=opt.twin{i}.tol.misGbRlx;
    end
    [G.Edges.combine,mGrains,twinGb] = ...
        MergeByBoundary(G,grains,tol,'Boundary Merged',opt);
    [G.Edges.combineRlx,mGrainsRlx,twinGbRlx] = ...
        MergeByBoundary(G,grains,tolRlx,'Boundary Merged Relaxed',opt);


    %Get the type of misorientation
    mori=inv(grains(G.Edges.pairs(:,1)).meanOrientation).*...
        grains(G.Edges.pairs(:,2)).meanOrientation; 
    tol=zeros(opt.nTwin,1);
    tolRlx=zeros(opt.nTwin,1);
    for i=1:opt.nTwin
        tol(i)=opt.twin{i}.tol.misMean;
        tolRlx(i)=opt.twin{i}.tol.misMeanRlx;
    end
    [~,G.Edges.meanType] = TestTwinRelationship(mori,tol,opt,G.Edges.type);
    [~,G.Edges.meanTypeRlx] = TestTwinRelationship(mori,tolRlx,opt,G.Edges.type);

    if opt.useMeanSmallGb
        [G.Edges.combine,G.Edges.GbInd] = ...
            UseMeanForSmallBoundary(G.Edges.pairs,G.Edges.GbInd,...
            G.Edges.combine,grains,G.Edges.meanType,opt);
    end
    
    if opt.mergeInclusion
        [G.Edges.type,G.Edges.combine] = ...
            AutoMergeInclusions(G.Edges.pairs,grains,G.Edges.combine,...
            opt.twinUnknown,G.Edges.meanTypeRlx);
    end
    
    %Overlayer graph on grains
    if opt.plot.do
        %Plot edge labeled graph
        labelNodes=false;labelEdges=opt.plot.labelEdges;plotG=true;legendOn=opt.plot.legendOn;
        fhandle = plotGraph(grains,[],G,...
            grains.meanOrientation,G.Nodes.Id,...
            labelNodes,labelEdges,legendOn,plotG,[]);
        mtexTitle('Full Graph')          
    end
end

function [type,combine] = AutoMergeInclusions(pairs,grains,combine,unknownType,meanTypeRlx)
    %Give inclusions an edge relationship regardless of whether they end up
    %being a twin. Usually they are twins or are 
    type=zeros(length(pairs),1,'int8');

    isInside = checkInside(grains, grains);
    [GrainIdInclusion,GrainIdWithInclusion] = find(isInside);
    
    epairs=[GrainIdInclusion,GrainIdWithInclusion];
    for i=1:length(GrainIdInclusion)
        eInd=all(epairs(i,:)==pairs,2) | all(fliplr(epairs(i,:))==pairs,2);
        if meanTypeRlx(i)~=0
            type(eInd)=meanTypeRlx(i);
        else
            type(eInd)=unknownType;
        end
        combine(eInd)=true;

    end
end


