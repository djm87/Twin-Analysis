function [G_clust,G] = FamilyVotes(G_clust,G,groupList,grains,mGrains,opt)
%FamilyVotes Computes a vote for twin/parent relationships 
%   Votes are made based on on the formula described in: 
%   Automatic twin statistcs from EBSD data 
%   Journal ofMicroscopy , Vol. 238, Pt 3 2010, pp. 218–229 
%   doi: 10.1111/j.1365-2818.2009.03343.x
%
%   One modification is performed, the effective schmid is per edge not per
%   family.
%
%   Inputs: the graph of clustered grains (G_clust) and weights (Sf, Af, and Bf)
%   Ouputs: Values needed to recalculate Votes are stored in G_clust along with 
%           along with votes using the weight passed in.
    if ~isempty(groupList)
        mineral=grains.mineral; %For defining phase boundary in the case there is unindexed data
        nArea=grains.area;
        %Loop over eachgroup of families
        for i=1:length(groupList)
            group=groupList(i);
            %Extract group to local variables for readability
            egroupId= find((group==G_clust.Edges.Group)); %converts logical arrays to indices
            if ~isempty(egroupId)
                ngroupId= find((group==G_clust.Nodes.Group));
                FID=G_clust.Nodes.FamilyID(ngroupId); %Family Id
                nID=G_clust.Nodes.Id(ngroupId); %Fragment grain Id
                nGb=grains(ngroupId); %Fragment grain boundary All time is spent here
        %         nGb=G_clust.Nodes.Gb(ngroupId); %Fragment grain boundary
                Area=nArea(ngroupId); %Fragment Area
        %         pairsGb=G_clust.Edges.Gb(egroupId); %Twin/parent pairs grain boundary
                pairs=G_clust.Edges.pairs(egroupId,1:2); %Twin/parent pairs
                eFamily=G_clust.Edges.FamilyID(egroupId,1:2);
                EffSF=G_clust.Edges.EffSF(egroupId,1:2); %Twin/parent Schmid factor
                type=G_clust.Edges.type(egroupId);
            
                %Family Areas to be stored in nodes 
                FArea= FamilyArea(Area,FID);
                G_clust.Nodes.FArea(ngroupId) = FArea;

                %Relative Family Areas to be stored in edges
                FRArea = AreaRatio(pairs,FArea,nID);
                G_clust.Edges.FRArea(egroupId,:) = FRArea; 

                %Family Boundaries (returns a cell)
                FgB = FamilyGrainBoundary(nGb,FID,nID,mineral);
    %             G_clust.Nodes.FgB(ngroupId) = FgB; %This is too large to store!

                %Boundary length ratio between connected families
                FRgB = GrainBoundaryRatio(pairs,FgB,eFamily,nID,mineral);
                G_clust.Edges.FRgB(egroupId,:) = FRgB;

                %Schmid factor difference (per edge not per family!) 
                FREffSF = SchmidFactorDifference(EffSF,pairs);  
                G_clust.Edges.FREffSF(egroupId,:) = FREffSF;

                %Calculate Vote (per edge)
                G_clust.Edges.Vote(egroupId,:) = CalcVote(FREffSF,FRArea,FRgB,type,opt);
            end
        end 

        %Transfer results to full graph
        G.Edges.FREffSF(:,:)=0;
        G.Edges.FREffSF(G.Edges.combineCleaned,:)=G_clust.Edges.FREffSF;
        G.Edges.FRgB(:,:)=0;
        G.Edges.FRgB(G.Edges.combineCleaned,:)=G_clust.Edges.FRgB;
        G.Edges.FRArea(:,:)=0;
        G.Edges.FRArea(G.Edges.combineCleaned,:)=G_clust.Edges.FRArea;
        G.Nodes.FArea(:)=0; %reset
        G.Nodes.FArea(G_clust.Nodes.Id)=G_clust.Nodes.FArea;
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

function [FgB] = FamilyGrainBoundary(nGb,FID,nID,mineral)
    gBId=nGb.boundary(mineral,mineral).grainId; %change this
    FgB=zeros(size(gBId,1),size(gBId,2));
    for j=1:length(nID)
        FgB(nID(j)==gBId)=FID(j);
    end
%         FgB=zeros(size(gBgId,1),size(gBgId,2),'logical');

%         FgB=cell(length(FID),1);
%         for j=1:max(FID)
% % %             nGbl=nGb(FID==j);
% % %             nGbF=nGbl{1};
% % %             for k=2:length(nGbl)
% % %                 nGbF=[nGbF;nGbl{k}];
% % %             end
%             nGbF=nGb(FID==j).boundary;
%             FgB(FID==j)={nGbF};
%         end   

end

function FRgB = GrainBoundaryRatio(pairs,FgB,eFamily,nID,mineral)
%     figure; plot(grains(nodeID),grains.meanOrientation(nodeID)); hold on;
%     gBColor=['k';'r';'g';'b';'y';'c'];
    for j=1:size(pairs,1)
%             n1=pairs(j,1);
%             n2=pairs(j,2);
        hasn1=any(eFamily(j,1)==FgB,2);
        hasn2=any(eFamily(j,2)==FgB,2);
        has0=any(0==FgB,2);
        n1gBLength=sum(hasn1);
        n2gBLength=sum(hasn2);
        n12gBLength=sum((~has0 & hasn1));
%             n1gB=FgB{n1==nID}(mineral,mineral);   
%             n2gB=FgB{n2==nID}(mineral,mineral);
%             n1gBLength=length(n1gB);   
%             n2gBLength=length(n2gB);%FgB{n2==nID}.length;
%             n1Ebsd=n1gB.ebsdId;
%             n2Ebsd=n2gB.ebsdId;
%             [~,loc]=intersect(n1Ebsd,n2Ebsd,'rows');
%             n12gBLength=length(loc);
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

function Vote = CalcVote(FREffSF,RFArea,FRgB,type,opt)
        nEdges=length(type);
        w=zeros(nEdges,3);
        for j=1:nEdges
            w(j,:)=opt.twin{1}.voteWeights;     
        end
        Vote(:,1) = w(:,1).*FREffSF(:,1)+w(:,2).*RFArea(:,1)+w(:,3).*FRgB(:,1);
        Vote(:,2) = w(:,1).*FREffSF(:,2)+w(:,2).*RFArea(:,2)+w(:,3).*FRgB(:,2);    
end

