function [grains] = TransferGtoGrains(G,grains)
%The graph structure is only need for computations, not for storing the
%analysis for datamining key variables which are outputed for further analysis
    grains.prop.Group = G.Nodes.Group;
    grains.prop.FamilyID = G.Nodes.FamilyID;
    grains.prop.EffSF = G.Nodes.EffSF;
    grains.prop.FArea = G.Nodes.FArea;
    grains.prop.typeUnknown = G.Nodes.typeUnknown;
    grains.prop.sigma13 = G.Nodes.sigma13;
    grains.prop.Generation = G.Nodes.Generation;
    grains.prop.ParentFamilyId = G.Nodes.ParentFamilyId;
    grains.prop.SchmidVariant = G.Nodes.SchmidVariant;
    grains.prop.K1NAng = G.Nodes.K1NAng;
    grains.prop.type = G.Nodes.type;
    grains.prop.twinThickness = G.Nodes.twinThickness;
end

