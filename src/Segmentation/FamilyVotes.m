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
    G_Family.Nodes.FgBL=zeros(length(G_Family.Nodes.Group),1);
    G_Family.Edges.FRArea=zeros(length(G_Family.Edges.Group),2);
    G_Family.Nodes.FVol=zeros(length(G_Family.Nodes.Group),1);
    G_Family.Edges.FRVol=zeros(length(G_Family.Edges.Group),2);
    G_Family.Edges.FRgB=zeros(length(G_Family.Edges.Group),2);
    G_Family.Edges.FSgB=zeros(length(G_Family.Edges.Group),1);
    G_Family.Edges.FREffSF=zeros(length(G_Family.Edges.Group),2);

    %Accessing the grain structure is very expensive and it it is better to
    %extract quantities in one go
    mineral=grains.mineral; %For defining phase boundary in the case there is unindexed data
    nArea=grains.area;
    if isfield(grains.prop,'volInitialODF')
        nVol=grains.prop.volInitialODF;
    else
        %Set so no influence in the vote
        nVol=ones(length(grains),1); 
    end
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
            Vol=nVol(ngroupId); %Grain volume in Initial ODF
            
            %Get Family quantities
            ngroupFId= find((group==G_Family.Nodes.Group)); %converts logical arrays to indices
            nFID=G_Family.Nodes.Family(ngroupFId);
            EffSF=G_Family.Edges.EffSF(egroupFId,1:2); %Twin/parent Schmid factor
            type=G_Family.Edges.meanType(egroupFId);
            FamilyID=G_Family.Edges.FamilyID(egroupFId,:);

            %Family Areas to be stored in nodes 
            [FArea,nFArea]= FamilyArea(Area,FID,nFID);
            G_Family.Nodes.FArea(ngroupFId) = nFArea;

            %Relative Family Areas to be stored in edges
            FRArea = AreaRatio(FamilyID,FArea);
            G_Family.Edges.FRArea(egroupFId,:) = FRArea; 

            %Max Family initial ODF volume
            
            [FVol,nFVol] = InitialODFVol(Vol,FID,nFID);
            G_Family.Nodes.FVol(ngroupFId) = nFVol;

            %Relative Family Initial ODF volume to be stored in edges
            FRVol = InitialODFVolRatio(FamilyID,FVol);
            G_Family.Edges.FRVol(egroupFId,:) = FRVol; 
            
            %Family Boundaries (returns a cell)
            [FgB] = FamilyGrainBoundary(nGrains,FID,nID,mineral);
            
            %Boundary length ratio between connected families
            [FRgB,FSgB,FgBL] = GrainBoundaryRatio(FgB,FamilyID);
            G_Family.Edges.FSgB(egroupFId) = FSgB;
            G_Family.Edges.FRgB(egroupFId,:) = FRgB;
            G_Family.Nodes.FgBL(ngroupFId) = FgBL; 
            
            %Schmid factor difference
            FREffSF = SchmidFactorDifference(EffSF,FamilyID);  
            G_Family.Edges.FREffSF(egroupFId,:) = FREffSF;

            %Calculate Vote (per edge)
            G_Family.Edges.Vote(egroupFId,:) = CalcVote(FREffSF,FRArea,FRgB,FRVol,type,opt);

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

function [FVol,nFVol] = InitialODFVol(Vol,FID,nFID)
    FVol=zeros(max(FID),1);
    nFVol=zeros(size(nFID,1),1);
    if all(Vol==1)
        FVol(:)=1;
        nFVol(:)=1;
    else
        for j=1:max(FID)
            tmp=max(Vol(FID==j));
            if ~isempty(tmp)
                FVol(j)=tmp;
                nFVol(nFID==j)=tmp;
            end
        end 
    end
end

function FRVol = InitialODFVolRatio(FamilyID,FVol)
    nFamilyPairs=size(FamilyID,1);
    FRVol=zeros(nFamilyPairs,2);
    for j=1:nFamilyPairs
        n1FVol=FVol(FamilyID(j,1));
        n2FVol=FVol(FamilyID(j,2));
        FRVol(j,1)=(n1FVol-n2FVol)/(n1FVol+n2FVol);
        FRVol(j,2)=(n2FVol-n1FVol)/(n1FVol+n2FVol);
    end 
end

function [FgB] = FamilyGrainBoundary(nGrains,FID,nID,mineral)
    gBId=nGrains.boundary(mineral,mineral).grainId; %change this
    FgB=zeros(size(gBId,1),size(gBId,2));
    for j=1:length(nID)
        hasgB=nID(j)==gBId;
        FgB(hasgB)=FID(j);
    end
end

function [FRgB,FSgB,FgBL] = GrainBoundaryRatio(FgB,FamilyID)
    nFamilyPairs=size(FamilyID,1);
    FRgB=zeros(nFamilyPairs,2);
    FSgB=zeros(nFamilyPairs,1);
    FgBL=zeros(max(max(FamilyID)),1);
    for j=1:nFamilyPairs
        hasn1=any(FamilyID(j,1)==FgB,2);
        hasn2=any(FamilyID(j,2)==FgB,2);
        has0=any(0==FgB,2);
        n1gBLength=sum(hasn1);
        n2gBLength=sum(hasn2);
        FgBL(FamilyID(j,1))=n1gBLength;
        FgBL(FamilyID(j,2))=n2gBLength;
        n12gBLength=sum((~has0 & hasn1));
        FSgB(j)=sum(hasn1&hasn2);
        FRgB(j,1)=(FSgB(j)/n2gBLength-FSgB(j)/n1gBLength);
        FRgB(j,2)=(FSgB(j)/n1gBLength-FSgB(j)/n2gBLength);
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

function Vote = CalcVote(FREffSF,RFArea,FRgB,FRVol,type,opt)
        nEdges=length(type);
        w=zeros(nEdges,length(opt.twin{1}.voteWeights));
        for j=1:nEdges
           if type(j)~=0 
              w(j,:)=opt.twin{type(j)}.voteWeights; 
           end
        end
        Vote(:,1) = w(:,1).*FREffSF(:,1)+w(:,2).*RFArea(:,1)+w(:,3).*FRgB(:,1)+w(:,4).*FRVol(:,1);
        Vote(:,2) = w(:,1).*FREffSF(:,2)+w(:,2).*RFArea(:,2)+w(:,3).*FRgB(:,2)+w(:,4).*FRVol(:,2);    
end

