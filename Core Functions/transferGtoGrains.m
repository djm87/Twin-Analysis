function [grains,mergedGrains] = transferGtoGrains(G,grains,mergedGrains,twin)
%The graph structure is only need for computations, not for storing the
%analysis for datamining key variables are outputed for further analysis
    grains.prop.EffSF=G.Nodes.EffSF;
    grains.prop.Group=G.Nodes.Group;
    grains.prop.FamilyID=G.Nodes.FamilyID;
    grains.prop.FArea=G.Nodes.FArea;
%     grains.prop.FgB=G.Nodes.FgB;
    grains.prop.Type=G.Nodes.Type;
    grains.prop.Generation=G.Nodes.Generation;
    grains.prop.twinCount=G.Nodes.twinCount;
    grains.prop.twinThickness=G.Nodes.twinThickness;

    grains.prop.schmidRank = cell(length(grains),1); %ranking of all twin variants
    grains.prop.schmidActive = zeros(length(grains),1); %active schmid value
    grains.prop.schmidActiveRank = zeros(length(grains),1); %active Rank
    grains.prop.schmidActiveN = zeros(length(grains),1); %Index
    grains.prop.schmid = cell(length(grains),1); %schmid value
    mergedGrains.prop.twinArea= zeros(length(mergedGrains),1); 
    
    area=grains.area;
    twinArea=zeros(length(area),1);
    twinArea(grains.prop.Type~=0)=area(grains.prop.Type~=0);
    openType=length(twin);
    groups=unique(G.Edges.Group);
    for i=1:length(groups)
        group=groups(i);
        %load edge and node properties for clustered fragments
        egroupId = find((group==G.Edges.Group)==true); %converts logical arrays to indices
        ngroupId = find((group==G.Nodes.Group)==true);
        typeKnown=~G.Nodes.typeUnknown(ngroupId);
        nId = G.Nodes.Id(ngroupId(typeKnown));
        nGeneration = -ones(length(nId),1,'int8');
        nType = -ones(length(nId),1,'int8');        
        eType = G.Edges.type(egroupId);
        etypeKnown=eType~=openType;
        eType=eType(etypeKnown);
        egroupId=egroupId(etypeKnown);
        eParent = G.Edges.Parent(egroupId,:);
        ePairs = G.Edges.pairs(egroupId,:);
        
        mergedGrains.prop.twinArea(group)=sum(twinArea(nId));
        grainId=ePairs(eParent);
        grains.prop.schmidRank(grainId) = G.Edges.schmidRank(egroupId);
        grains.prop.schmidActive(grainId) = G.Edges.schmidActive(egroupId);
        grains.prop.schmidActiveRank(grainId) = G.Edges.schmidActiveRank(egroupId);
        grains.prop.schmid(grainId) = G.Edges.schmid(egroupId);
        grains.prop.schmidActiveN(grainId) = G.Edges.schmidActiveN(egroupId);
    end


end

