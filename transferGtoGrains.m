function [grains] = transferGtoGrains(G_Complete,grains)
%The graph structure is only need for computations, not for storing the
%analysis for datamining key variables are outputed for further analysis
    grains.prop.EffSF=G_Complete.Nodes.EffSF;
    grains.prop.Group=G_Complete.Nodes.Group;
    grains.prop.FamilyID=G_Complete.Nodes.FamilyID;
    grains.prop.FArea=G_Complete.Nodes.FArea;
    grains.prop.FgB=G_Complete.Nodes.FgB;
    grains.prop.Type=G_Complete.Nodes.Type;
    grains.prop.Generation=G_Complete.Nodes.Generation;
    grains.prop.twinCount=G_Complete.Nodes.twinCount;
    grains.prop.twinThickness=G_Complete.Nodes.twinThickness;

    grains.prop.schmidRank = zeros(length(grains),1);
    grains.prop.schmidActive = zeros(length(grains),1);
    grains.prop.schmidActiveRank = zeros(length(grains),1);
    grains.prop.schmidActiveN = zeros(length(grains),1);
    grains.prop.schmid = zeros(length(grains),size(G_Complete.Edges.schmid,2));


    for i=1:max(G_Complete.Edges.Group) 
        egroupId = find((i==G_Complete.Edges.Group)==true); %converts logical arrays to indices
        ePairs = G_Complete.Edges.pairs(egroupId,:);
        eParent = G_Complete.Edges.Parent(egroupId,:);
        
        grains.prop.schmidRank(ePairs(eParent)) = G_Complete.Edges.schmidRank(egroupId);
        grains.prop.schmidActive(ePairs(eParent)) = G_Complete.Edges.schmidActive(egroupId);
        grains.prop.schmidActiveRank(ePairs(eParent)) = G_Complete.Edges.schmidActiveRank(egroupId);
        grains.prop.schmid(ePairs(eParent),:) = G_Complete.Edges.schmid(egroupId,:);
        grains.prop.schmid(ePairs(eParent)) = G_Complete.Edges.schmidActiveN(egroupId);

    end

end

