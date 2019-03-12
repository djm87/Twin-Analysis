function [twinVF,totalArea,totalTwinVF] = getTwinFractions(G_Complete,grains,twin)
%getTwinFractions loops through the available modes and gives area
%fractions and the total area
    %Extract the grains areas
    area = G_Complete.Nodes.Area;
    totalArea = sum(area);
    
    %Initialize twin fractions
    ntwin = length(twin);   
    twinVF = zeros(ntwin,1);
    
    %compute twin fractions
    for i=1:ntwin
        twinVF(i) = sum(area(G_Complete.Nodes.Type == i)) / totalArea;
    end
    
    totalTwinVF =sum(twinVF) / totalArea;
    
    
end

