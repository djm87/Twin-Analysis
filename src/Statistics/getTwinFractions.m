function [mGrains,twinVF] = getTwinFractions(G,grains,mGrains,opt)
%getTwinFractions loops through the available modes and gives area
%fractions and the total area

    %Extract the grain areas
    area=grains.area;
    mArea=mGrains.area;
    twinVF.totalArea = sum(area);
    
    %Initialize twin fractions  
    nGen=max(G.Nodes.Generation);
    twinVF.perMode = zeros(opt.nTwin,1);
    twinVF.perModeAndGen = zeros(opt.nTwin,nGen);
    
    %compute combined statistics
    for i=1:opt.nTwin
        for j=1:nGen
            twinVF.perModeAndGen(i,j)=...
                sum(area(G.Nodes.type == i & G.Nodes.Generation == j)) / twinVF.totalArea;
        end
        twinVF.perMode(i) = sum(twinVF.perModeAndGen(i,:));
    end
    twinVF.total = sum(twinVF.perMode);
    
    %compute the per merged grain twin information 
    TVF_perGrainPerModeAndGen = zeros(length(mGrains),opt.nTwin,nGen);
    TVF_perGrainPerMode = zeros(length(mGrains),opt.nTwin);
    for i=1:length(mGrains)
        nId=find(G.Nodes.Group==i);
        nType=G.Nodes.type(nId);
        nGeneration=G.Nodes.Generation(nId);
        nArea=area(nId);
        nMArea=mArea(i);
        
        for j=1:opt.nTwin
            for k=1:nGen
                TVF_perGrainPerModeAndGen(i,j,k)=...
                    sum(nArea(nType == j & nGeneration == k)) / nMArea;
            end
            TVF_perGrainPerMode(i,j)=sum(TVF_perGrainPerModeAndGen(i,j,:));
        end
    end
    
    mGrains.prop.TVF_perGrainPerModeAndGen=TVF_perGrainPerModeAndGen;
    mGrains.prop.TVF_perGrainPerMode=TVF_perGrainPerMode;


end

