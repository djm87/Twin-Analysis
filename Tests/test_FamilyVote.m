
% for i=1:max(G_clustClean.Edges.Group)
%     egroupId= find((i==G_clustClean.Edges.Group)==true); %convert logical arrays to indices
%     ngroupId= find((i==G_clustClean.Nodes.Group)==true);
% 
%     ori=G_clustClean.Nodes.meanOrientation(ngroupId);
%     G_clustClean.Nodes.FamilyID(ngroupId)=GetFamily(ori)
% end


Vote1={};
Vote2={};
G_clustClean.Nodes.FamilyVote=zeros(length(G_clustClean.Nodes.Id),1);
for i=1:max(G_clustClean.Edges.Group)
    egroupId= find((i==G_clustClean.Edges.Group)==true); %convert logical arrays to indices
    ngroupId= find((i==G_clustClean.Nodes.Group)==true);
    Area=G_clustClean.Nodes.Area(ngroupId)
    FamilyID=G_clustClean.Nodes.FamilyID(ngroupId)
    nodeID=G_clustClean.Nodes.Id(ngroupId)
    nodeGb=G_clustClean.Nodes.Gb(ngroupId)
    Gb_pairs=G_clustClean.Edges.Gb(egroupId)
    pairs=G_clustClean.Edges.pairs(egroupId,1:2)
    SF=G_clustClean.Edges.SF(egroupId,1:2)
    
    %Family Areas to be stored in nodes 
    FA=zeros(length(FamilyID),1);
    for j=1:max(FamilyID)
        FA(FamilyID==j)=sum(Area(FamilyID==j))
    end
    
    %Relative Family Areas to be stored in edges
    clear 'RFA1' 'RFA2'
    for j=1:size(pairs,1)
        n1=pairs(j,1)
        n2=pairs(j,2)
        n1AF=FA(n1==nodeID)
        n2AF=FA(n2==nodeID)
        RFA1(j)=(n1AF-n2AF)/(n1AF+n2AF)
        RFA2(j)=(n2AF-n1AF)/(n1AF+n2AF)
    end
    
    %Family Boundaries
    FgB=cell(length(FamilyID),1);
    for j=1:max(FamilyID)
        nGb=nodeGb(FamilyID==j)
        
        nGbF=nGb{1};
        for k=2:length(nGb)
            nGbF=[nGbF;nGb{k}];
        end
        FgB(FamilyID==j)={nGbF};
    end  
    
    %Boundary length ratio between connected families
%     figure; plot(grains(nodeID),grains.meanOrientation(nodeID)); hold on;
%     gBColor=['k';'r';'g';'b';'y';'c'];
    clear 'FRgB1' 'FRgB2'
    for j=1:size(pairs,1)
        n1=pairs(j,1);
        n2=pairs(j,2);
        n1gB=FgB{n1==nodeID}('Mg','Mg');   
        n2gB=FgB{n2==nodeID}('Mg','Mg');
        n1gBLength=FgB{n1==nodeID}.length;   
        n2gBLength=FgB{n2==nodeID}.length;
        n1Ebsd=n1gB.ebsdId;
        n2Ebsd=n2gB.ebsdId;
        [boundaryEbsdId,loc]=intersect(n1Ebsd,n2Ebsd,'rows');
        n12gBLength=n1gB(loc).length;
%         plot(n1gB(loc),gBColor(j),'linecolor',gBColor(j),'linewidth',3)
        FRgB1(j)=(n12gBLength/n2gBLength-n12gBLength/n1gBLength);
        FRgB2(j)=(n12gBLength/n1gBLength-n12gBLength/n2gBLength);
    end
%     hold off
    
    %Schmid factor difference 
    clear 'FRSF1' 'FRSF2'
    for j=1:size(pairs,1)
        FRSF2(j)=diff(SF(j,1:2),1); %diff defined as second-first
        FRSF1(j)=-FRSF2(j);
    end    
    Sf=1; Af=1; Bf=1;
    Vote1{i}=Sf*FRSF1+Af*RFA1+Bf*FRgB1
    Vote2{i}=Sf*FRSF2+Af*RFA2+Bf*FRgB2
end
%%
G_clustClean.Nodes.Properties.UserData.Mineral=grains.mineral
w=[1,1,1]
G_wihtVoteTest = FamilyVotes(G_clustClean,w)
%%
i=3
indx=3
Vote1{indx}
Vote2{indx}
G_wihtVoteTest.Edges.Vote(egroupId,:)