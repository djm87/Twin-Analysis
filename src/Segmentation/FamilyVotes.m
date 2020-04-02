function [G_Family] = FamilyVotes(G_Family,G_clust,groupList,grains,opt)
%FamilyVotes Computes a vote for twin/parent relationships 
%   Votes are made based on on the formula described in: 
%   Automatic twin statistcs from EBSD data 
%   Journal ofMicroscopy , Vol. 238, Pt 3 2010, pp. 218–229 
%   doi: 10.1111/j.1365-2818.2009.03343.x
%
%   Inputs: the graph of Families and weights for each vote (per twin type)
%   Ouputs: Values needed to recalculate Votes are stored in G_Family along with 
%           along with votes using the weights passed in.

    %If there is nothing to do return
    if isempty(groupList)
       return 
    end

    %Initialize quantities that are being computed and stored
    G_Family.Nodes.FArea=zeros(length(G_Family.Nodes.Group),1);
    G_Family.Edges.FRArea=zeros(length(G_Family.Edges.Group),2);
    G_Family.Edges.FRgB=zeros(length(G_Family.Edges.Group),2);
    G_Family.Edges.FREffSF=zeros(length(G_Family.Edges.Group),2);
    
    %Accessing the grain structure is very expensive and it it is better to
    %extract quantities in one go
    mineral=grains.mineral; %For defining phase boundary in the case there is unindexed data
    nArea=grains.area;
    
    %Loop over eachgroup of families
    for i=1:length(groupList)
        
        group=groupList(i);
        %Extract group to local variables for readability
        egroupFId= find((group==G_Family.Edges.Group));
        if ~isempty(egroupFId)
            %Get cluster quantities
            ngroupId= find((group==G_clust.Nodes.Group));
            FID=G_clust.Nodes.FamilyID(ngroupId); %Family Id
            nID=G_clust.Nodes.Id(ngroupId); %Fragment grain Id
            
            %Extract a subset of grains and areas. For performance
            nGrains=grains(ngroupId); %Fragment grain boundary
            Area=nArea(ngroupId); %Fragment Area
            
            %Get Family quantities
            ngroupFId= find((group==G_Family.Nodes.Group)); %converts logical arrays to indices
            nFID=G_Family.Nodes.Family(ngroupFId);
            EffSF=G_Family.Edges.EffSF(egroupFId,1:2); %Twin/parent Schmid factor
            type=G_Family.Edges.meanTypeRlx(egroupFId);
            FamilyID=G_Family.Edges.FamilyID(egroupFId,:);

            %Family Areas to be stored in nodes 
            [FArea,nFArea]= FamilyArea(Area,FID,nFID);
            G_Family.Nodes.FArea(ngroupFId) = nFArea;

            %Relative Family Areas to be stored in edges
            FRArea = AreaRatio(FamilyID,FArea);
            G_Family.Edges.FRArea(egroupFId,:) = FRArea; 

            %Family Boundaries (returns a cell)
            FgB = FamilyGrainBoundary(nGrains,FID,nID,mineral);
%             G_clust.Nodes.FgB(ngroupId) = FgB; %This is too large to store!

            %Boundary length ratio between connected families
            FRgB = GrainBoundaryRatio(FgB,FamilyID);
            G_Family.Edges.FRgB(egroupFId,:) = FRgB;

            %Schmid factor difference (per edge not per family!) 
            FREffSF = SchmidFactorDifference(EffSF,FamilyID);  
            G_Family.Edges.FREffSF(egroupFId,:) = FREffSF;

            %Calculate Vote (per edge)
            G_Family.Edges.Vote(egroupFId,:) = CalcVote(FREffSF,FRArea,FRgB,type,opt);
        end
    end 
end

function [FArea,nFArea] = FamilyArea(Area,FID,nFID)
    FArea=zeros(max(FID),1);
    nFArea=zeros(size(nFID,1),1);
    for j=1:max(FID)
        FArea(j)=sum(Area(FID==j));
        nFArea(nFID==j)=FArea(j);
    end 
end

function RFArea = AreaRatio(FamilyID,FArea)
    nFamilyPairs=size(FamilyID,1);
    RFArea=zeros(nFamilyPairs,2);
    for j=1:nFamilyPairs
        n1FArea=FArea(FamilyID(j,1));
        n2FArea=FArea(FamilyID(j,2));
        RFArea(j,1)=(n1FArea-n2FArea)/(n1FArea+n2FArea);
        RFArea(j,2)=(n2FArea-n1FArea)/(n1FArea+n2FArea);
    end 
end

function [FgB] = FamilyGrainBoundary(nGrains,FID,nID,mineral)
    gBId=nGrains.boundary(mineral,mineral).grainId; %change this
    FgB=zeros(size(gBId,1),size(gBId,2));
    for j=1:length(nID)
        FgB(nID(j)==gBId)=FID(j);
    end
end

function FRgB = GrainBoundaryRatio(FgB,FamilyID)
    nFamilyPairs=size(FamilyID,1);
    FRgB=zeros(nFamilyPairs,2);
    for j=1:nFamilyPairs
        hasn1=any(FamilyID(j,1)==FgB,2);
        hasn2=any(FamilyID(j,2)==FgB,2);
        has0=any(0==FgB,2);
        n1gBLength=sum(hasn1);
        n2gBLength=sum(hasn2);
        n12gBLength=sum((~has0 & hasn1));
        FRgB(j,1)=(n12gBLength/n2gBLength-n12gBLength/n1gBLength);
        FRgB(j,2)=(n12gBLength/n1gBLength-n12gBLength/n2gBLength);
    end

end

function FREffSF= SchmidFactorDifference(EffSF,FamilyID)
    %Note if unknown type then EffSF will be 0 and no contribution to
    %the vote will occure for that edge.
    nFamilyPairs=size(FamilyID,1);
    FREffSF=zeros(nFamilyPairs,2);
    for j=1:nFamilyPairs
        FREffSF(j,2)=diff(EffSF(j,1:2)); %diff defined as second-first
        FREffSF(j,1)=-FREffSF(j,2);
    end     
end

function Vote = CalcVote(FREffSF,RFArea,FRgB,type,opt)
        nEdges=length(type);
        w=zeros(nEdges,3);
        for j=1:nEdges
            w(j,:)=opt.twin{1}.voteWeights;     
        end
        Vote(:,1) = w(:,1).*FREffSF(:,1)+w(:,2).*RFArea(:,1)+w(:,3).*FRgB(:,1);
        Vote(:,2) = w(:,1).*FREffSF(:,2)+w(:,2).*RFArea(:,2)+w(:,3).*FRgB(:,2);    
end

