function fhandle = plotGraph(grains,mergedGrains,graph,grainPlot,nId,labelNodes,labelEdges,legendOn,plotG,options)
%plotGraph plots a graph overlay on grains 
%Inputs: grains (full 2dgrain object from EBSD
%        mergedGrains (contains the grouped grain boundary corresponding to grains)
%        graph (the paired down graph object that is to be plotted)
%        grainPlot (the scalar or orientaitons values to be be plotted (length nId)
%        nId (the node ids of grains to be ploted)
    fhandle=figure;
    plot(grains(nId),...
        grainPlot(nId),'Micronbar','off','silent');
    hold on 
    if plotG
        p=plot(graph,'XData',graph.Nodes.centroids(:,1),...
                'YData',graph.Nodes.centroids(:,2),'displayName','graph');
    end
    if isempty(options)
        p.EdgeColor='k';p.MarkerSize=2;p.Marker='s';p.NodeColor='k';
        mGlc='k'; mGlw=2;p.NodeFontSize=18;p.EdgeFontSize=14;p.LineWidth=1
    else
        p.EdgeColor=options{1};p.MarkerSize=options{2};p.Marker=options{3};p.NodeColor=options{4};p.EdgeFontSize=options{5};
        mGlc=options{6}; mGlw=options{7};
    end
    if ~isempty(mergedGrains)
        plot(mergedGrains.boundary,'linecolor',mGlc,'linewidth',mGlw,...
            'linestyle','-','displayName','merged grains')    
    end
    
    if labelEdges && plotG
        pairs1=graph.Edges.pairs(:,1);
        pairs2=graph.Edges.pairs(:,2);
        for i=1:length(nId)
            pairs1(pairs1==graph.Nodes.Id(i))=i;
            pairs2(pairs2==graph.Nodes.Id(i))=i;
        end
        labeledge(p,pairs1,...
            pairs2,graph.Edges.GlobalID);
    end
    
    if ~labelNodes
        p.NodeLabel={}
    end
    
    if ~legendOn
        legend off
    end
    hold off
    
end

