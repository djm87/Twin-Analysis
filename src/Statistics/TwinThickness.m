function [G] = TwinThickness(G,grains,opt)
%TwinThickness Takes the ellipsoids fit to twin grains and corrects it by
%the cosine of the angle between k1 axis and out of plane direction 
    [~,~,b_min] = fitEllipse(grains);
    twinThickness = zeros(length(grains),1);
    for i =1:length(grains)
        nId = G.Nodes.Id(i);
        nType = G.Nodes.type(i);
        if nType~=opt.twinUnknown && nType>0
            if length(opt.twin{nType}.k1)==1
                twinThickness(nId)=b_min(nId)*abs(cos(G.Nodes.K1NAng(nId)));
            else
                twinThickness(nId)=b_min(nId);
            end
        end
    end
    G.Nodes.twinThickness = twinThickness;
    
end