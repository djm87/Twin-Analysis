function twinBoundary = GetTwinBoundaries(G_Complete,grains,twin,Mistol)
%Returns the twin bounaries
    ntwins=length(twin);
    mineral=grains.mineral;
    gB=grains.boundary;
    gB_mineral = gB(mineral,mineral);
    twinBoundary={};
    for i=1:ntwins    
        isTwinning = angle(gB_mineral.misorientation,twin{i}.RMT) < Mistol;
        twinBoundary{i} = gB_mineral(isTwinning);
    end
    figure; 
    plot(grains,G_Complete.Nodes.Type,'Micronbar','off')
    hold on 

    colors=hsv(ntwins);
    for i=1:ntwins
        plot(twinBoundary{i},'linecolor',colors(i,:),'linewidth',2,'DisplayName',twin{i}.name);
    end
    hold off
end

