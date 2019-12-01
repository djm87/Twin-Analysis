function [G,mergedGrains,twinBoundary] = MergeByBoundary(G,grains,Mistol,twin)
%Applies standard mtex grain merging from boundaries.

    ntwins=length(twin);
    mineral=grains.mineral;
    gB=grains.boundary;
    gB_mineral = gB(mineral,mineral);
    twinBoundary={};
    cnt=0;
    for i=1:ntwins 
        twinBoundary{i}=[];
        for j=1:length(twin{i}.variantsToUse) %for double twins
            cnt=cnt+1;
            isTwinning = angle(gB_mineral.misorientation,twin{i}.RMT(j)) < Mistol;
            twinBoundary{cnt} = gB_mineral(isTwinning);
        end
    end
    combinedTwinBoundary=vertcat(twinBoundary{:});
    [mergedGrains,parentId] = merge(grains,combinedTwinBoundary);
    
%     %Reduce edges to clusters that will be tested using mean relaxed
%     for i=1:length(mergedGrains) 
%         grainid=find(parentId == mergedGrains(i).id)
%         if length(grainid)>1
%             
%         end
%     end
    
    pairs=G.Edges.pairs;
    combineBoundary=cell(length(mergedGrains),1);
    parfor i=1:length(mergedGrains) 
%         grainid = grains(parentId == mergedGrains(i).id).id;
        grainid=find(parentId == mergedGrains(i).id);
        if length(grainid)>1
            combine1=zeros(length(pairs(:,1)),1,'logical');
            combine2=zeros(length(pairs(:,1)),1,'logical');
            
            %for each merged grain, identify all edges that have the grains
            for j=1:length(grainid)
                combine1(find(grainid(j)==pairs(:,1)))=true;
                combine2(find(grainid(j)==pairs(:,2)))=true;
            end
            %Sum the pairs to see which edges are in the merged
            %grains. This creates a combine that is based on boundary
            combineBoundary{i}=find((combine1+combine2)==2);
        end
    end
    G.Edges.combineBoundary=zeros(length(G.Edges.pairs),1,'logical');
    G.Edges.combineBoundary(vertcat(combineBoundary{:}))=true;
    %Map the merged grains to pairs
%     G.Edges.combineMerge=zeros(length(G.Edges.pairs),1);
%     G.Edges.groupMerge=zeros(length(G.Edges.pairs),1);
%     for i=1:length(G.Edges.pairs)
%           nodes=G.Edges.pairs(i,:);
%           if all(ismember(nodes,grains(combine).id))
%              G.Edges.combineMerge(i)=true;
%              G.Edges.groupMerge(i)=group(nodes(1));
%           end
%     end
       
%     figure; 
%     plot(grains,grains.meanOrientation,'Micronbar','off')
%     hold on 
%     colors=hsv(ntwins);
%     for i=1:ntwins
%         plot(twinBoundary{i},'linecolor',colors(i,:),'linewidth',2,'displayName',twin{i}.name);
%     end
%     plot(mergedGrains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-',...
%         'displayName','merged grains')
%     p=plot(G,'XData',G.Nodes.centroids(:,1),...
%             'YData',G.Nodes.centroids(:,2),'displayName','graph');
%         hold off
%         p.Marker='s';p.NodeColor='k';p.MarkerSize=3;p.EdgeColor='k';
%     hold off
    
    

end

