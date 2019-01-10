function twinBoundary = GetTwinBoundaries(G_Complete,grains,twins,twinNames,Mistol)
%Returns the twin bounaries
    ntwins=length(twins);
    mineral=grains.mineral;
    gB=grains.boundary;
    gB_mineral = gB(mineral,mineral);
    twinBoundary={};
    for i=1:ntwins    
        isTwinning = angle(gB_mineral.misorientation,twins{i}) < Mistol;
        twinBoundary{i} = gB_mineral(isTwinning);
    end
    figure; 
    plot(grains(G_Complete.Nodes.Id(G_Complete.Nodes.Type>-2)),...
        G_Complete.Nodes.Type(G_Complete.Nodes.Type>-2),'Micronbar','off')
    hold on 
    colors=hsv(ntwins);
    for i=1:ntwins
        plot(twinBoundary{i},'linecolor',colors(i,:),'linewidth',2,'DisplayName',twinNames{i});
    end
    hold off
end

