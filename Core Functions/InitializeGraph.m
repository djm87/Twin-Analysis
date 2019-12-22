function [G,time] = InitializeGraph(ebsd,grains,twin,Mistol,meanMistol,...
    meanMistolRelaxed,doPlot,doLabel,time)
    %This function initializes graph objects for grains. 
    %Initialization entails computation of all local grain properties such as
    %orientation, centroid, aspectRatio, etc.. and is stored in edge pair
    %(intergranular) or nodes (grains). 

    tic

    %Compute grain neighbors
    [counts,pairs] = neighbors(grains);

    %Initialize graph
    s=pairs(:,1);
    t=pairs(:,2);
    G=graph(s,t);

    % Add node labels and other grain level information 
    G.Nodes.Name=cellstr(int2str(grains.id));
    G.Nodes.Id=str2num(cell2mat(G.Nodes.Name));
    G.Nodes.Area=grains.area; 
    G.Nodes.Perimeter=grains.perimeter; 
    G.Nodes.AspectRatio=grains.aspectRatio;
    G.Nodes.Paris=grains.paris;
    G.Nodes.centroids=grains.centroid;
    G.Nodes.meanOrientation=grains.meanOrientation;
    G.Nodes.Properties.UserData.mineral=grains.mineral; %For single phase material

    % Add intergranular information
    G.Edges.pairs=pairs; %this contains the grain ids corresponding to grains.id
    G.Edges.GlobalID=[1:size(pairs,1)]';
    G.Edges.Parent=zeros(size(pairs,1),2);
    
    %label internal grains twin type in last twin (i.e. type=length(twin))
    type=zeros(length(G.Edges.pairs),1,'int8');
    G.Nodes.typeUnknown=zeros(length(G.Nodes.Id),1,'logical');

    isInside = checkInside(grains, grains);
    [GrainIdInclusion,GrainIdWithInclusion] = find(isInside);
    G.Nodes.typeUnknown(GrainIdInclusion)=true;
    epairs=[GrainIdInclusion,GrainIdWithInclusion];
    for i=1:length(GrainIdInclusion)
        eId=find(all(epairs(i,:)==G.Edges.pairs,2) | all(fliplr(epairs(i,:))==G.Edges.pairs,2));
        type(eId)=length(twin);
    end

    %Test if mean misorientation is a twin type so we can cluster grains
%     [combine,type] = TestTwinRelationship(mori,meanMistol,twin,type);

    G.Edges.type=type; %Twin relation type (# from twin list definitions)
%     G.Edges.combine=combine; %whether pairs should be grouped into grains

    %Overlayer graph on grains
    if doPlot==true
        %Plot edge labeled graph
        labelNodes=false;labelEdges=doLabel;plotG=true;legendOn=false;
        fhandle = plotGraph(grains,[],G,...
            grains.meanOrientation,G.Nodes.Id,...
            labelNodes,labelEdges,legendOn,plotG,[]);
        mtexTitle('Initial Graph')          
    end

    if ~isfield(time,'InitializeGraph')
        time.InitializeGraph=0;
    end
    time.InitializeGraph=time.InitializeGraph+toc;
end

