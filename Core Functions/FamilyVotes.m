function G = FamilyVotes(G,grains,w)
%FamilyVotes Computes a vote for twin/parent relationships 
%   Votes are made based on on the formula described in: 
%   Automatic twin statistcs from EBSD data 
%   Journal ofMicroscopy , Vol. 238, Pt 3 2010, pp. 218–229 
%   doi: 10.1111/j.1365-2818.2009.03343.x
%
%   One modification is performed, the effective schmid is per edge not per
%   family.
%
%   Inputs: the graph of clustered grains (G) and weights (Sf, Af, and Bf)
%   Ouputs: Values needed to recalculate Votes are stored in G along with 
%           along with votes using the weight passed in.

    %Initialize arrays in G 
    G.Nodes.FArea=zeros(length(G.Nodes.Area),1);
    G.Nodes.FgB=cell(length(G.Nodes.Area),1);
    
    G.Edges.FRArea=zeros(length(G.Edges.EffSF),2);
    G.Edges.FRgB=zeros(length(G.Edges.EffSF),2);
    G.Edges.FREffSF=zeros(length(G.Edges.EffSF),2);
    G.Edges.Vote=zeros(length(G.Edges.EffSF),2);
    
    mineral=G.Nodes.Properties.UserData.mineral; %For defining phase boundary in the case there is unindexed data
    
    %Loop over eachgroup of families
    for i=1:max(G.Edges.Group)
        %Extract group to local variables for readability
        egroupId= find((i==G.Edges.Group)==true); %converts logical arrays to indices
        ngroupId= find((i==G.Nodes.Group)==true);
        FID=G.Nodes.FamilyID(ngroupId); %Family Id
        nID=G.Nodes.Id(ngroupId); %Fragment grain Id
        nGb=grains(ngroupId); %Fragment grain boundary
%         nGb=G.Nodes.Gb(ngroupId); %Fragment grain boundary
        Area=G.Nodes.Area(ngroupId); %Fragment Area
%         pairsGb=G.Edges.Gb(egroupId); %Twin/parent pairs grain boundary
        pairs=G.Edges.pairs(egroupId,1:2); %Twin/parent pairs
        EffSF=G.Edges.EffSF(egroupId,1:2); %Twin/parent Schmid factor
        if ~isempty(pairs)
            %Family Areas to be stored in nodes 
            FArea= FamilyArea(Area,FID);
            G.Nodes.FArea(ngroupId) = FArea;

            %Relative Family Areas to be stored in edges
            FRArea = AreaRatio(pairs,FArea,nID);
            G.Edges.FRArea(egroupId,:) = FRArea; 

            
            %Family Boundaries (returns a cell)
            FgB = FamilyGrainBoundary(nGb,FID);
            G.Nodes.FgB(ngroupId) = FgB;

            %Boundary length ratio between connected families
            FRgB = GrainBoundaryRatio(FgB,pairs,nID,mineral);
            G.Edges.FRgB(egroupId,:) = FRgB;

            %Schmid factor difference (per edge not per family!) 
            FREffSF = SchmidFactorDifference(EffSF,pairs);  
            G.Edges.FREffSF(egroupId,:) = FREffSF;

            %Calculate Vote (per edge)
            G.Edges.Vote(egroupId,:) = CalcVote(FREffSF,FRArea,FRgB,w);
        end
    end    
    
    function FArea = FamilyArea(Area,FID)
        FArea=zeros(length(FID),1);
        for j=1:max(FID)
            FArea(FID==j)=sum(Area(FID==j));
        end 
    end

    function RFArea = AreaRatio(pairs,FArea,nID)
        for j=1:size(pairs,1)
            n1=pairs(j,1);
            n2=pairs(j,2);
            n1FArea=FArea(n1==nID);
            n2FArea=FArea(n2==nID);
            RFArea(j,1)=(n1FArea-n2FArea)/(n1FArea+n2FArea);
            RFArea(j,2)=(n2FArea-n1FArea)/(n1FArea+n2FArea);
        end 
    end

    function FgB = FamilyGrainBoundary(nGb,FID)
        FgB=cell(length(FID),1);
        for j=1:max(FID)
%             nGbl=nGb(FID==j);
%             nGbF=nGbl{1};
%             for k=2:length(nGbl)
%                 nGbF=[nGbF;nGbl{k}];
%             end
            nGbF=nGb(FID==j).boundary;
            FgB(FID==j)={nGbF};
        end   
    end

    function FRgB = GrainBoundaryRatio(FgB,pairs,nID,mineral)
    %     figure; plot(grains(nodeID),grains.meanOrientation(nodeID)); hold on;
    %     gBColor=['k';'r';'g';'b';'y';'c'];
        for j=1:size(pairs,1)
            n1=pairs(j,1);
            n2=pairs(j,2);
            n1gB=FgB{n1==nID}(mineral,mineral);   
            n2gB=FgB{n2==nID}(mineral,mineral);
            n1gBLength=FgB{n1==nID}.length;   
            n2gBLength=FgB{n2==nID}.length;
            n1Ebsd=n1gB.ebsdId;
            n2Ebsd=n2gB.ebsdId;
            [boundaryEbsdId,loc]=intersect(n1Ebsd,n2Ebsd,'rows');
            n12gBLength=n1gB(loc).length;
    %         plot(n1gB(loc),gBColor(j),'linecolor',gBColor(j),'linewidth',3)
            FRgB(j,1)=(n12gBLength/n2gBLength-n12gBLength/n1gBLength);
            FRgB(j,2)=(n12gBLength/n1gBLength-n12gBLength/n2gBLength);
        end
    %     hold off
    end

    function FREffSF= SchmidFactorDifference(EffSF,pairs)
        %Note if unknown type then EffSF will be 0 and no contribution to
        %the vote will occure for that edge.
        for j=1:size(pairs,1)
            FREffSF(j,2)=diff(EffSF(j,1:2)); %diff defined as second-first
            FREffSF(j,1)=-FREffSF(j,2);
        end     
    end

    function Vote = CalcVote(FREffSF,RFArea,FRgB,w)
            Vote(:,1) = w(1)*FREffSF(:,1)+w(2)*RFArea(:,1)+w(3)*FRgB(:,1);
            Vote(:,2) = w(1)*FREffSF(:,2)+w(2)*RFArea(:,2)+w(3)*FRgB(:,2);    
    end
end

