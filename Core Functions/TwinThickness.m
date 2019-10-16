function [G_Complete] = TwinThickness(G_Complete,grains,twin)
%TwinThickness Takes the ellipsoids fit to twin grains and corrects it by
%the cosine of the angle between k1 axis and out of plane direction 
    [~,~,b_min] = fitEllipse(grains);
    twinThickness = zeros(length(grains),1);
    for i =1:length(grains)
        nId = G_Complete.Nodes.Id(i);
        nType = G_Complete.Nodes.Type(i);
        if nType>0
            if length(twin{nType}.k1)==1
                twinThickness(nId)=b_min(nId)*G_Complete.Nodes.k1NormalAngle(nId);
            else
                twinThickness(nId)=b_min(nId);
            end
        end
    end
    G_Complete.Nodes.twinThickness = twinThickness;
    
end

