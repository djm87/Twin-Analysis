function [G_Family,iouterGroup,exitExternal] = GraphEditor(groups,iouterGroup,G_Family,G_clust,G,grains,mergedGrains,value,plotNeighbors,plotEdgeLabel,plotMergedGrainId,plotGrainId,mode)
%ClusterFilter loops through clusters in prefiltered groups, allowing a 
%user to add/remove nodes edges and relationships
%Inputs:  groups - list of merged grain Ids 
%         G - graph
%         grains - 2dgrain object
%         value - matches length grains and contains a scalar value or orientation
%Outputs: Most outputs go to .txt files used in ClusterGrains.
%         The Graph object is updated if an edge is added.
%         runCleanupAgain is used in CleanFamilyTree to recluster grains
%         ind is returned for whether CleanFamilyTree should move on to the
%         next cluster.

    exitExternal=false;
    igroup=1
    while igroup<length(groups)+1
        group=groups(igroup);
        G_clust_sub = subgraph(G_clust,find(group==G_clust.Nodes.Group));

        % Mode 1 is the family graph editor
        if mode==1
            G_Family.Edges.eGlobalID=[1:numedges(G_Family)]';
            G_Family_sub = subgraph(G_Family,find(group==G_Family.Nodes.Group));
            C = centrality(G_Family_sub,'hubs');
            
            %Plot cluster
            h=newMtexFigure;
            plotFamilyId=true;
            plotGrainId=false;
            plotEdgeLabel=false;
            plotClusterGraph(group,G_clust_sub,G_clust,...
                grains,value,mergedGrains,plotGrainId,plotNeighbors,...
                plotMergedGrainId,plotEdgeLabel,plotFamilyId);            

            nextAxis;
            plotFamilyGraph(full(adjacency(G_Family_sub)));
%             mtexColorbar('Location','eastoutside','Ticks',1:3)
            colorbar(h.children(1),'Location','eastoutside')  
            set(h.parent,'WindowStyle','docked')
              
            %Call the Family editor
            [G_Family_sub,G_Family,igroup,iouterGroup,exitExternal,...
                mode,value,plotNeighbors]=familyEditorInterface(group,igroup,iouterGroup,...
                G_Family,G_Family_sub,G_clust,grains,exitExternal,mode,value); 
            
            %Cleanup the figure
            close(h.parent);
            
            plotFamilyId=false;
            plotGrainId=true;
            plotEdgeLabel=true;
        %Mode 2 is the clust graph editor
        elseif mode ==2
            %Plot cluster
            h=figure;
            plotFamilyId=false;
            plotClusterGraph(group,G_clust_sub,G_clust,...
                grains,value,mergedGrains,plotGrainId,plotNeighbors,...
                plotMergedGrainId,plotEdgeLabel,plotFamilyId);
            set(h,'WindowStyle','docked')
            
            %Call the cluster editor
            [igroup,iouterGroup,exitExternal,plotNeighbors,value,mode,plotEdgeLabel,...
                plotGrainId,plotMergedGrainId] = clusterEditorInterface(group,...
                igroup,iouterGroup,G_clust_sub,G_clust,G,grains,exitExternal,...
                plotNeighbors,value,mode,plotEdgeLabel,plotGrainId,plotMergedGrainId); 
            
            %Cleanup the figure
            close(h);
            
        end
        if exitExternal
           return;
        end
        
    end
    fprintf('Finished in GraphEditor\n')
end
function plotFamilyGraph(FamilyMatrix)
    G_Family_clust=digraph(FamilyMatrix);
%     set(0,'DefaultFigureWindowStyle','docked')
%     set(0,'DefaultFigureWindowStyle','docked');
%     warning('off', 'MATLAB:Figure:SetPosition');
    p=plot(G_Family_clust);
%     layout(p,'force','Iterations',30)
    p.EdgeFontSize=13;p.NodeFontSize=14;p.ArrowSize=15;p.EdgeColor=[0,0,0];
    p.NodeColor=[1,0,0];p.MarkerSize=8;
    pair1=G_Family_clust.Edges.EndNodes(:,1);
    pair2=G_Family_clust.Edges.EndNodes(:,2);
    npairs=length(pair1);
    labeledge(p,pair1,pair2,1:npairs); 
%     set(0,'DefaultFigureWindowStyle','normal');
%     warning('off', 'MATLAB:uitools:uimode:callbackerror');
%     warning('on', 'MATLAB:Figure:SetPosition');
%     G_Family_clust=flipedge(G_Family_clust,3)
%     G_Family_clust = rmedge(G_Family_clust,3)
%     FamilyMatrix = full(adjacency(G_Family_clust))
end
function [G_Family_sub,G_Family,igroup,iouterGroup,exitExternal,...
                mode,value,plotNeighbors]=familyEditorInterface(group,igroup,iouterGroup,...
                G_Family,G_Family_sub,G_clust,grains,exitExternal,mode,value) 
    try
    fprintf('Edge List \n')
    for j=1:numedges(G_Family_sub)
        fprintf('Id: %5d, type: %3d, Pair/Parent: %5d %5d,  %1d %1d\n',j,...
            G_Family_sub.Edges.meanType(j),G_Family_sub.Edges.pairs(j,:),G_Family_sub.Edges.Parent(j,:))
    end
    
    %Give options to perform
    fprintf('===================================\n')
    fprintf('processing options\n')
    fprintf('Press..\n')
    fprintf('enter to proceed to next grain cluster\n')
    fprintf('0 to exit editor and return exit flag\n')
    fprintf('1 to remove an edge\n')
    fprintf('2 to remove edges connected to Family\n')
    fprintf('3 to flip parent relationship\n')
    fprintf('4 to plot the grain a different way\n')
    fprintf('5 to switch editor mode to cluster\n')
    fprintf('6 go back to previous grain\n')
    fprintf('7 cycle\n')
    
    %Initialize
    plotNeighbors=false;
    
    %Get main operation to perform
    while true
        inputMsg='Enter Number: ';
        errorMsg='Please input a valid option or hit enter to go to the next grain';
        option = getUserInput(inputMsg,errorMsg,'scalar',true);
        
        if option >8 | option < 0
            fprintf('Enter a valid option\n')
        elseif option==0
            fprintf('Exiting GraphEditor\n')
            exitExternal=true;
            return
        elseif option==1 
            while true
                inputMsg='Enter edges to remove: ';
                errorMsg='Must specify a scalar 1 or a vector [1,2] or hit enter to return to menu';
                eID_sub = getUserInput(inputMsg,errorMsg,'vector',true);
                if ~any(eID_sub<1)
                    eID=G_Family_sub.Edges.eGlobalID(eID_sub)
                    G_Family=rmedge(G_Family,eID);
                    return;
                elseif isempty(nodeId)
                    break;
                else
                    fprintf('%s\n',errorMsg)
                end
            end
        elseif option==2 
            while true
                inputMsg='Enter families to remove edges from: ';
                errorMsg='Must specify a scalar 1 or a vector [1,2] or hit enter to return to menu';
                nID_sub = getUserInput(inputMsg,errorMsg,'vector',true);
                if ~any(nID_sub<1)
                    toRemove=zeros(numedges(G_Family_sub),1,'logical');
                    for i=1:length(nID_sub)
                        toRemove(any(G_Family_sub.Edges.FamilyID==nID_sub(i),2))=true;
                    end
                    eID=G_Family_sub.Edges.eGlobalID(toRemove)
                    G_Family=rmedge(G_Family,eID);
                    return;
                elseif isempty(nodeId)
                    break;
                else
                    fprintf('%s\n',errorMsg)
                end
            end
        elseif option==3 
            while true
                inputMsg='Enter edges to flip: ';
                errorMsg='Must specify a scalar 1 or a vector [1,2] or hit enter to return to menu';
                eID_sub = getUserInput(inputMsg,errorMsg,'vector',true);
                if ~any(eID_sub<1)
                    eID=G_Family_sub.Edges.eGlobalID(eID_sub)
                    G_Family=flipedge(G_Family,eID);
                    return;
                elseif isempty(nodeId)
                    break;
                else
                    fprintf('%s\n',errorMsg)
                end
            end
%         elseif option==4
%             while true
%                 inputMsg='Enter edge to add: ';
%                 errorMsg='Must specify a scalar 1 or a vector [1,2] or hit enter to return to menu';
%                 ePairs = getUserInput(inputMsg,errorMsg,'vector',true);
%                 if ~any(eID_sub<1) && all(size(ePairs)==[1,2])
%                     G_Family=addedge(G_Family,ePairs(1),ePairs(2));
%                     return;
%                 elseif isempty(nodeId)
%                     break;
%                 else
%                     fprintf('%s\n',errorMsg)
%                 end
%             end
        elseif option==4
            fprintf('1 to plot mean orientation\n')
            fprintf('2 to plot FamilyId\n')
            fprintf('3 to plot EffSchmid\n')
            fprintf('4 to plot Initial ODF VF\n')
            while true
                inputMsg='Enter plot option number: ';
                errorMsg='Please input a valid option or hit enter to go back to menu';
                plotOption = getUserInput(inputMsg,errorMsg,'scalar',true);
                %Check the type of input
                if isempty(plotOption)
                    break;
                elseif (plotOption==1 || plotOption==2 || plotOption==3 || plotOption==4)
                    if plotOption==1
                        value=grains.meanOrientation;
                    elseif plotOption==2
                        value=G_clust.Nodes.FamilyID;
                    elseif plotOption==3
                        inputMsg='Enter twin type id: ';
                        errorMsg='Please input a valid twin type id';
                        twinType = getUserInput(inputMsg,errorMsg,'scalar',false);
                        value=G_clust.Nodes.EffSF(:,twinType);
                    elseif plotOption==4
                        try
                            value=grains.prop.volInitialODF;
                        catch
                            fprintf('It seems the volume fraction from the initial ODF is not available\n') 
                        end
                    end
                    inputMsg='Enter if neighbor should be plotted (0 or 1): ';
                    errorMsg='Please input a valid option';
                    plotNeighbors= getUserInput(inputMsg,errorMsg,'logical',false);
                    break;
                else 
                    fprintf('Please enter a valid plot option\n')
                end
            end
            exitExternal=false;
            return
        elseif option==5
            mode=2;
            return;
        elseif option==6
            fprintf('going back to previous grain\n');
            igroup=igroup-1;
            iouterGroup=iouterGroup-1;
            if igroup==0
                fprintf('At first grain, cannot go back\n');  
                igroup=1;
            else
                return
            end
        elseif option==7
            fprintf('Cycle GraphEditor\n');
            return
        elseif isempty(option)
            fprintf('Moving on to the next grain\n');
            igroup=igroup+1;
            iouterGroup=iouterGroup+1;
            return
        else
            %This should never happen.. but if it does catch it
            error('Error: unhandled grain')
        end
    end
    catch
       fprintf('Your intput caused an error.. no worries your changes are safe') 
    end
end
function plotClusterGraph(group,G_clust_sub,G_clust,...
    grains,value,mergedGrains,plotGrainId,plotNeighbors,...
    plotMergedGrainId,plotEdgeLabel,plotFamilyId)
    
    nId_sub = G_clust_sub.Nodes.Id;
    eGlobalId_sub = G_clust_sub.Edges.GlobalID;
    
    %if plot neighbors 
    mGrainInd=group;
    if plotNeighbors
        [~,pairsTmp]=neighbors(mergedGrains);
        isPair=pairsTmp(:,1)==group | pairsTmp(:,2)==group;
        allGroupsInd=unique(pairsTmp([isPair,isPair]));
        nId=intersect_wrepeat(allGroupsInd,G_clust.Nodes.Group);
    else
       nId=nId_sub;
       allGroupsInd=mGrainInd;
    end

    %Plot grain cluster
%     set(0,'DefaultFigureWindowStyle','docked');
%     warning('off', 'MATLAB:Figure:SetPosition');
    plot(grains(nId),value(nId),'noBoundary')
%         mtexColorMap(hsv)
    
    if plotFamilyId
        text(grains(nId),int2str(G_clust.Nodes.FamilyID(nId)));
    end

    if plotGrainId
        text(grains(nId),int2str(nId));
    end

    if plotMergedGrainId
        grainName=cell(length(allGroupsInd),1);
        for ii=1:length(allGroupsInd)
            grainName{ii}=sprintf('M%d',mergedGrains(allGroupsInd(ii)).id);
        end
        text(mergedGrains(allGroupsInd),grainName(:));        
    end

    hold on
    if ~isempty(mergedGrains)
        if plotNeighbors
            plot(mergedGrains(allGroupsInd).boundary,'linecolor','k','linewidth',2,...
                'linestyle','-','displayName','merged grains')
        end
        plot(mergedGrains(mGrainInd).boundary,'linecolor','w','linewidth',3,...
            'linestyle','-','displayName','merged grains')
    end

    if plotEdgeLabel
%         toremove=ones(length(G_clust.Nodes.Id),1,'logical');
%         toremove(unique(G_clust.Edges.pairs))=false;
% %         toremove(nId)=false;
%         G_Removed=rmnode(G_clust,find(toremove));
% 
%         [ii,jj]=unique(G_Removed.Edges.GlobalID,'stable');
%         n=numel(eGlobalId);
%         pos = zeros(n,1);
%         for k = 1:n
%              pos(k)=jj(find(ii == eGlobalId(k)));
%         end
%         toremove=ones(size(G_Removed.Edges.GlobalID,1),1,'logical');
%         toremove(pos)=false;
%         G_Removed=rmedge(G_Removed,find(toremove));

        p=plot(G_clust_sub,'XData',G_clust_sub.Nodes.centroids(:,1),...
            'YData',G_clust_sub.Nodes.centroids(:,2),'displayName','graph');
        p.EdgeFontSize=13;
        p.EdgeColor=[0,0,0];
        p.LineWidth=1;
        pairs1=G_clust_sub.Edges.EndNodes(:,1);
        pairs2=G_clust_sub.Edges.EndNodes(:,2);
%         for j=1:length(G_Removed.Nodes.Id)
%             pairs1(pairs1==G_clust_sub.Nodes.Id(j))=j;
%             pairs2(pairs2==G_clust_sub.Nodes.Id(j))=j;
%         end
        if ~isempty(pairs1)
            labeledge(p,pairs1,...
                pairs2,G_clust_sub.Edges.GlobalID); 
        end
        nNodes=length(unique(G_clust_sub.Nodes.Id));
        labelnode(p,1:nNodes,strings(nNodes,1));
    end
    hold off   
    
%     set(0,'DefaultFigureWindowStyle','normal');
%     warning('off', 'MATLAB:uitools:uimode:callbackerror');
%     warning('on', 'MATLAB:Figure:SetPosition');
end
function [igroup,iouterGroup,exitExternal,plotNeighbors,value,mode,plotEdgeLabel,...
    plotGrainId,plotMergedGrainId] = clusterEditorInterface(group,...
    igroup,iouterGroup,G_clust_sub,G_clust,G,grains,exitExternal,...
    plotNeighbors,value,mode,plotEdgeLabel,plotGrainId,plotMergedGrainId) 
    
    try
    nFamily=G_clust_sub.Nodes.FamilyID;
    nId_sub=G_clust_sub.Nodes.Id;
    
    %give node and edge info
    fprintf('===================================\n')
    fprintf('Group %d\n',group)
    fprintf('Node List \n')
    for j=1:max(nFamily)
        nId_family=nId_sub(j==nFamily);
        fprintf('Family %d, Node Id ',j)
        for k=1:length(nId_family)
            fprintf('%d ',nId_family(k))
        end
        fprintf('\n')
    end
    fprintf('Edge List \n')

    for j=1:numedges(G_clust_sub)
        fprintf('Id: %5d, type: %3d Node, Pair/Parent: %5d %5d,  %1d %1d\n',G_clust_sub.Edges.GlobalID(j),...
            G_clust_sub.Edges.type(j),G_clust_sub.Edges.pairs(j,:),G_clust_sub.Edges.Parent(j,:))
    end

    %Give options to perform
    fprintf('===================================\n')
    fprintf('processing options\n')
    fprintf('Press..\n')
    fprintf('enter to proceed to next grain cluster\n')
    fprintf('0 to exit editor and return exit flag\n')
    fprintf('1 to remove an edge\n')
    fprintf('2 to add an edge\n')
    fprintf('3 to remove all edges connected to a node\n')
    fprintf('4 to try adding all edges connected to a node\n')
    fprintf('5 merge merged grains\n')
    fprintf('6 to plot the grain a different way\n')
    fprintf('7 to change labeling\n')
    fprintf('8 to get misorientation of grains\n')
    fprintf('9 to get all edge relationships with node\n')
    fprintf('10 go back to previous grain\n')
    fprintf('11 cycle\n')
    fprintf('12 switch to Family editor\n')
    
    %Get main operation to perform

    while true
        inputMsg='Enter Number: ';
        errorMsg='Please input a valid option or hit enter to go to the next grain';
        option = getUserInput(inputMsg,errorMsg,'scalar',true);
        
        if option >12 | option < 0
            fprintf('Enter a valid option\n')
        elseif option==0
            fprintf('exiting GraphEditor\n')
            exitExternal=true;
            return
        elseif option == 1
            while true
                inputMsg='Enter edge Id: ';
                errorMsg='Please input a valid edge id or hit enter to return to menu';
                edgeId = getUserInput(inputMsg,errorMsg,'vector',true);
                if ~any(edgeId==0)
                   break;
                elseif isempty(edgeId)
                   break;
                else
                    fprintf('%s\n',errorMsg)
                end
            end

            fid = fopen('eRemoveList.txt', 'a+');
            for j=1:length(edgeId)
                fprintf(fid, '%d\n', edgeId(j));
            end
            fclose(fid);
            
        elseif option == 2                    
            while true
                inputMsg='Enter node pair for edge: ';
                errorMsg='Must specify a vector of two nodes [1,2] or for multiple [12,1;1,4] or hit enter to return to menu';
                ePair = getUserInput(inputMsg,errorMsg,'vector',true);
                if all(~any(ePair==0,2))
                    break;
                elseif isempty(ePair)
                    break;
                else
                    fprintf('%s\n',errorMsg)
                end
            end

            fid = fopen('eAddList.txt', 'a+');
            for j=1:size(ePair,1)
                fprintf(fid, '%d %d\n', ePair(j,1), ePair(j,2));
            end
            fclose(fid);
        elseif option == 3                   
            while true
                inputMsg='Enter list of node ids';
                errorMsg='Must specify a sclar 1 or a vector [1,2] or hit enter to return to menu';
                nodeId = getUserInput(inputMsg,errorMsg,'vector',true);
                if ~any(nodeId==0)
                    break;
                elseif isempty(nodeId)
                    break;
                else
                    fprintf('%s\n',errorMsg)
                end
            end
            fid = fopen('nRemoveList.txt', 'a+');
            for j=1:length(nodeId)
                fprintf(fid, '%d\n', nodeId(j));
            end
            fclose(fid);
        elseif option == 4                   
            while true
                inputMsg='Enter list of node ids';
                errorMsg='Must specify a sclar 1 or a vector [1,2] or hit enter to return to menu';
                nodeId = getUserInput(inputMsg,errorMsg,'vector',true);
                if ~any(nodeId==0)
                    break;
                elseif isempty(nodeId)
                    break;
                else
                    fprintf('%s\n',errorMsg)
                end
            end
            fid = fopen('nAddList.txt', 'a+');
            for j=1:length(nodeId)
                fprintf(fid, '%d\n', nodeId(j));
            end
            fclose(fid);
        elseif option == 5  
            while true
                inputMsg='Enter merge grain pair for merging: ';
                errorMsg='Must specify a vector of two grains [1,2] or for multiple [12,1;1,4] or hit enter to return to menu';
                mgPairs = getUserInput(inputMsg,errorMsg,'vector',true);
                if all(~any(mgPairs==0,2))
                    break;
                elseif isempty(mgPairs)
                    break;
                else
                    fprintf('%s\n',errorMsg)
                end
            end

            %Find the pairs and add to eAddList.txt
            ePairs=[];
            for j=1:size(mgPairs,1)
                ngroupId1 = find((mgPairs(j,1)==G_clust.Nodes.Group)==true);
                ngroupId2 = find((mgPairs(j,2)==G_clust.Nodes.Group)==true);
                gBId2=grains(ngroupId2).boundary.grainId;
                for k=1:length(ngroupId1)
                    [row,~]=find(ngroupId1(k)==gBId2);
                    uniqueNeighbors=unique(gBId2(row,:));
                    uniqueNeighbors(uniqueNeighbors==ngroupId1(k))=[];
                    ePairs=vertcat(ePairs,[uniqueNeighbors,ngroupId1(k)*ones(length(uniqueNeighbors),1)]);
                end
            end     
            fid = fopen('eAddList.txt', 'a+');
            for j=1:size(ePairs,1)
                fprintf(fid, '%d %d\n', ePairs(j,1), ePairs(j,2));
            end
            fclose(fid);

        elseif option == 6
            fprintf('1 to plot mean orientation\n')
            fprintf('2 to plot FamilyId\n')
            fprintf('3 to plot EffSchmid\n')
            while true
                inputMsg='Enter plot option number: ';
                errorMsg='Please input a valid option or hit enter to go back to menu';
                plotOption = getUserInput(inputMsg,errorMsg,'scalar',true);
                %Check the type of input
                if isempty(plotOption)
                    break;
                elseif (plotOption==1 || plotOption==2 || plotOption==3)
                    if plotOption==1
                        value=grains.meanOrientation;
                    elseif plotOption==2
                        value=G_clust.Nodes.FamilyID;
                    elseif plotOption==3
                        inputMsg='Enter twin type id: ';
                        errorMsg='Please input a valid twin type id';
                        twinType = getUserInput(inputMsg,errorMsg,'scalar',false);
                        try
                            value=G_clust.Nodes.EffSF(:,twinType);
                        catch
                            fprintf('Invalid twin type specified\n')
                            break;
                        end
                    end
                    inputMsg='Enter if neighbor should be plotted (0 or 1): ';
                    errorMsg='Please input a valid option';
                    plotNeighbors= getUserInput(inputMsg,errorMsg,'logical',false);
                    break;
                else 
                    fprintf('Please enter a valid plot option\n')
                end
            end
            exitExternal=false;
            return
        elseif option == 7
            fprintf('1 plot Neighbors\n')
            fprintf('2 plot EdgeLabel\n')
            fprintf('3 plot GrainId\n')
            fprintf('4 plot Merged GrainId\n')

            while true
                inputMsg='Enter plot option number: ';
                errorMsg='Please input a valid option or hit enter to return to menu';
                plotOption = getUserInput(inputMsg,errorMsg,'scalar',true);
                %Check the type of input
                if isempty(plotOption)
                    break;
                elseif (plotOption==1 || plotOption==2 || plotOption==3 || plotOption==4)
                    if plotOption==1
                        inputMsg='On or off (0 or 1): ';
                        errorMsg='Please input a valid option';
                        plotNeighbors = getUserInput(inputMsg,errorMsg,'logical',true);
                    elseif plotOption==2
                        inputMsg='On or off (0 or 1): ';
                        errorMsg='Please input a valid option';
                        plotEdgeLabel = getUserInput(inputMsg,errorMsg,'logical',true);
                    elseif plotOption==3
                        plotGrainId=true;
                        plotMergedGrainId=false;
                    elseif plotOption==4                        
                        plotGrainId=false;
                        plotMergedGrainId=true;                          
                    end
                end
            end
            return
        elseif option==8   
            while true
                inputMsg='Enter node pair for comparing: ';
                errorMsg='Must specify a vector of two grains [1,2] or for multiple [12,1;1,4] or hit enter to return to menu';
                compList = getUserInput(inputMsg,errorMsg,'vector',true);
                if all(~any(compList==0,2))
                    break;
                elseif isempty(compList)
                    break;
                else
                    fprintf('%s\n',errorMsg)
                end
            end

            angle(grains(compList(:,1)).meanOrientation, grains(compList(:,2)).meanOrientation)./ degree
        elseif option==9 
            while true
                inputMsg='Enter nodes to get edge relationships for: ';
                errorMsg='Must specify a sclar 1 or a vector [1,2] or hit enter to return to menu';
                nodeId = getUserInput(inputMsg,errorMsg,'vector',true);
                if ~any(nodeId==0)
                    break;
                elseif isempty(nodeId)
                    break;
                else
                    fprintf('%s\n',errorMsg)
                end
            end

            for i=1:length(nodeId)
               edgeId=find(G.Edges.pairs(:,1)==nodeId(i) | G.Edges.pairs(:,2)==nodeId(i)); 
               mori=angle(grains(G.Edges.pairs(edgeId,1)).meanOrientation,...
                   grains(G.Edges.pairs(edgeId,2)).meanOrientation)./ degree;
               fprintf('Relationships for Node %d\n',nodeId(i))
               for j=1:length(edgeId)
                fprintf('EdgeID: %6d, NodeIds: (%5d,%5d), Gb Rlx Twin Type: %d, Misorientation: %6.2f\n',...
                    edgeId(j),G.Edges.pairs(edgeId(j),1),G.Edges.pairs(edgeId(j),2),G.Edges.GBTypeRlx(edgeId(j)),mori(j)) 
               end
               fprintf('\n')
            end
            
        elseif option==10
            fprintf('going back to previous grain\n');
            igroup=igroup-1;
            iouterGroup=iouterGroup-1;
            if igroup==0
                fprintf('At first grain, cannot go back\n');  
                igroup=1;
            else
                return
            end
        elseif option==11
            fprintf('Cycle GraphEditor\n');
            return
        elseif option==12
            mode=1;
            return
        elseif isempty(option)
            fprintf('Moving on to the next grain\n');
            igroup=igroup+1;
            iouterGroup=iouterGroup+1;
            return
        else
            %This should never happen.. but if it does catch it
            error('Error: unhandled grain')
        end
    end
    catch
       fprintf('Your intput caused an error.. no worries your changes are safe') 
    end
end
function val=getUserInput(inputMsg,errorMsg,type,canBeEmpty)
    try
    switch type
        case {'scalar','Scalar'}   
            while true
                try
                    val=input(inputMsg);
                catch Error
                    disp(Error.message)
                    fprintf('%s\n',errorMsg)
                end
                if isnumeric(val)
                    break;
                elseif isempty(val)
                    if canBeEmpty
                        break;
                    else
                       fprintf('Input is not allowed to be empty\n') 
                    end
                end
                fprintf('%s\n',errorMsg)
            end
        case {'vector','Vector'}   
            while true
                try
                    val=input(inputMsg);
                catch Error
                    disp(Error.message)
                    fprintf('%s\n',errorMsg)
                end
                if isvector(val) || ismatrix(val)
                    break
                elseif isempty(val)
                    if canBeEmpty
                        break;
                    else
                       fprintf('Input is not allowed to be empty\n') 
                    end
                end
            end
        case {'logical','Logical'}   
            while true
                try
                    val=input(inputMsg);
                    val=logical(val);
                catch Error
                    disp(Error.message)
                    fprintf('%s\n',errorMsg)
                end
                if islogical(val)
                    break;
                elseif isempty(val)
                    if canBeEmpty
                        break;
                    else
                       fprintf('Input is not allowed to be empty\n') 
                    end
                end
                fprintf('%s\n',errorMsg)
            end
    end
    catch
       fprintf('Your intput caused an error.. no worries your changes are safe') 
    end
end