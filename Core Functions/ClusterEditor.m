function [G,runCleanupAgain,ind] = ClusterEditor(groups,G,grains,mergedGrains,value,ind,plotNeighbors,enforceClusterOnlyMod)
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

    rEdge=[];
    runCleanupAgain=false;
    for i=1:length(groups) 
        group=groups(i);
        egroupId = find((group==G.Edges.Group)==true); %converts logical arrays to indices
        ngroupId = find((group==G.Nodes.Group)==true);
        nId = G.Nodes.Id(ngroupId);
        nFamily = G.Nodes.FamilyID(ngroupId);
        eGlobalId = G.Edges.GlobalID(egroupId);  
        
        %if plot neighbors 
        if plotNeighbors
            [~,pairsTmp]=neighbors(grains(nId));
            nId=unique(pairsTmp);
            allGroups=unique(G.Nodes.Group(nId));
            nId=[];
            for j=1:length(allGroups)
                ngroupIdAll = find((allGroups(j)==G.Nodes.Group)==true);
                nId = vertcat(nId,G.Nodes.Id(ngroupIdAll));
            end
        end
        
        %Plot grain cluster
        set(0,'DefaultFigureWindowStyle','docked')
        h=figure; plot(grains(nId),value(nId))
        text(grains(nId),int2str(nId))
        hold on
        if ~isempty(mergedGrains)
            plot(mergedGrains(group).boundary,'linecolor','k','linewidth',2,...
                'linestyle','-','displayName','merged grains')
        end
        toremove=ones(length(G.Nodes.Id),1,'logical');
        toremove(unique(G.Edges.pairs(egroupId,:)))=false;
%         toremove(nId)=false;
        G_Removed=rmnode(G,find(toremove));
        p=plot(G_Removed,'XData',G_Removed.Nodes.centroids(:,1),...
            'YData',G_Removed.Nodes.centroids(:,2),'displayName','graph');
        pairs1=G_Removed.Edges.pairs(:,1);
        pairs2=G_Removed.Edges.pairs(:,2);
        for j=1:length(G_Removed.Nodes.Id)
            pairs1(pairs1==G_Removed.Nodes.Id(j))=j;
            pairs2(pairs2==G_Removed.Nodes.Id(j))=j;
        end
        if ~isempty(pairs1)
            labeledge(p,pairs1,...
                pairs2,G_Removed.Edges.GlobalID); 
        end
        hold off   
        set(0,'DefaultFigureWindowStyle','normal')
        %give node and edge info
        fprintf('===================================\n')
        fprintf('Group %d\n',group)
        fprintf('Node List \n')
        for j=1:max(nFamily)
            nId_family=nId(j==nFamily);
            fprintf('Family %d, Node Id ',j)
            for k=1:length(nId_family)
                fprintf('%d ',nId_family(k))
            end
            fprintf('\n')
        end
        fprintf('Edge List \n')
        for j=1:length(egroupId)
            fprintf('Id: %5d, Node Pair: %5d %5d\n',G.Edges.GlobalID(egroupId(j)),G.Edges.pairs(egroupId(j),:))
        end

        %Give options to perform
        fprintf('===================================\n')
        fprintf('processing options\n')
        fprintf('Press..\n')
        fprintf('enter to proceed to next grain cluster\n')
        fprintf('0 to exit with editor\n')
        fprintf('1 to specify node as twin and parent of none\n')
        fprintf('2 to specify node as child of none\n')
        fprintf('3 to set an edge relationship\n')
        fprintf('4 to remove an edge\n')
        fprintf('5 to add an edge\n')
        fprintf('6 to remove all edges connected to a node\n')
        fprintf('7 to try adding all edges connected to a node\n')
        %Get main operation to perform
        nextGrain=false;
        while true
            while true
                option=input('enter number: ');
                if isempty(option)
                    break;
                elseif length(option)>1
                    %cycle
                elseif isnumeric(option) && (option >= 0 && option < 5 )
                    break;
                end
                fprintf('Please input a valid option\n')
            end
            if option==0
                fprintf('exiting ClusterEditor\n')
                close(h);
                return;
            elseif option==1
                while true
                    try
                        fprintf('Specify a number or vector e.g. 1 or [1,2]\n')
                        nodeId=input('enter list of Node Id: ');
                    catch Error
                        disp(Error.message)
                        fprintf('Must specify a number or vector e.g. 1 or [1,2]\n')
                    end
                    %Check the type of input
                    if ~any(nodeId==0) && isvector(nodeId) 
                        %Check that the input makes sense
                        flagBreak=true;
                        for j=1:length(nodeId) 
                           if and(~any(nodeId(j)==nId),enforceClusterOnlyMod)
                              flagBreak=false; 
                              fprintf('You specified a nodeId outside the cluster\n')
                           end
                        end
                        if flagBreak
                            break;
                        end
                    elseif nodeId==0 
                        break;
                    end
                end
                fid = fopen('notParent.txt', 'a+');
                for j=1:length(nodeId)
                    fprintf(fid, '%d\n', nodeId(j));
                end
                fclose(fid);
            elseif option == 2
                while true
                    try
                        nodeId=input('enter list of Node Id: ');
                    catch Error
                        disp(Error.message)
                        fprintf('Must specify a number or vector e.g. 1 or [1,2]\n')
                    end
                    %Check the type of input
                    if ~any(nodeId==0) && isvector(nodeId)
                        %Check that the input makes sense
                        flagBreak=true;
                        for j=1:length(nodeId) 
                           if and(~any(nodeId(j)==nId),enforceClusterOnlyMod)
                              flagBreak=false; 
                              fprintf('You specified a nodeId outside the cluster\n')
                           end
                        end
                        if flagBreak
                            break;
                        end
                    end
                end
                fid = fopen('notTwin.txt', 'a+');
                for j=1:length(nodeId)
                    fprintf(fid, '%d\n', nodeId(j));
                end
                fclose(fid);
            elseif option == 3  
                fprintf('To specify a relation enter edge id, and logical pair e.g. [124 1 0] where 1 is the parent\n')
                while true
                    try
                        eRelation=input('specify edge relation: ');
                    catch Error
                        disp(Error.message)
                        fprintf('Must specify a vector e.g. [124 1 0]\n')
                    end
                    if eRelation(1)~=0 && isvector(eRelation)
                        %Check that the input makes sense
                        flagBreak=true;

                       if and(~any(eRelation(1)==eGlobalId),enforceClusterOnlyMod)
                          flagBreak=false; 
                          fprintf('You specified an edge outside the cluster\n')
                       end

                       if sum(eRelation(2:3))~=1
                          flagBreak=false; 
                          fprintf('You need to specified a parent twin relationship as 1 0 or 0 1 where 1 is the parent\n')
                       end
                        if flagBreak
                            break;
                        end
                    end
                end
                fid = fopen('eRelation.txt', 'a+');
                fprintf(fid, '%d %d %d\n', eRelation);
                fclose(fid);
                
            elseif option == 4                    
                while true
                    try
                        edgeId=input('enter list of edge Id: ');
                    catch Error
                        disp(Error.message)
                        fprintf('Must specify a number or vector e.g. 1 or [1,2]\n')
                    end
                    if ~any(edgeId==0) && isvector(edgeId)
                        %Check that the input makes sense
                        flagBreak=true;
                        for j=1:length(edgeId) 
                           if and(~any(edgeId(j)==eGlobalId), enforceClusterOnlyMod)
                              flagBreak=false; 
                              fprintf('You specified a edge outside the cluster\n')
                           end
                        end
                        if flagBreak
                            break;
                        end
                    end
                end

                fid = fopen('eRemoveList.txt', 'a+');
                for j=1:length(edgeId)
                    fprintf(fid, '%d\n', edgeId(j));
                end
                fclose(fid);
                rEdge=zeros(length(edgeId),1);
                for j=1:length(edgeId)
                    rEdge(j)=find(edgeId(j)==eGlobalId);
                end
                G=removeEdge(G,rEdge,egroupId);
                ind=ind+1;
                runCleanupAgain=true;
            elseif option == 5                    
                while true
                    try
                        ePair=input('enter node pair for edge: ');
                    catch Error
                        disp(Error.message)
                        fprintf('Must specify a vector of two nodes [1,2]\n')
                    end
                    if all(~any(ePair==0,2)) && (isvector(ePair) || ismatrix(ePair))
                        break;
                    else
                        fprintf('Specify a node id pairs e.g.[12,1] or for multiple [12,1;1,4]\n')
                    end
                end

                fid = fopen('eAddList.txt', 'a+');
                for j=1:size(ePair,1)
                    fprintf(fid, '%d %d\n', ePair(j,1), ePair(j,2));
                end
                fclose(fid);
            elseif option == 6                   
                while true
                    try
                        nodeId=input('enter list of Node Id: ');
                    catch Error
                        disp(Error.message)
                        fprintf('Must specify a number or vector e.g. 1 or [1,2]\n')
                    end
                    %Check the type of input
                    if ~any(nodeId==0) && isvector(nodeId)
                        %Check that the input makes sense
                        flagBreak=true;
                        for j=1:length(nodeId) 
                           if and(~any(nodeId(j)==nId),enforceClusterOnlyMod)
                              flagBreak=false; 
                              fprintf('You specified a nodeId outside the cluster\n')
                           end
                        end
                        if flagBreak
                            break;
                        end
                    end
                end
                fid = fopen('nRemoveList.txt', 'a+');
                for j=1:length(nodeId)
                    fprintf(fid, '%d\n', nodeId(j));
                end
                fclose(fid);
            elseif option == 7                   
                while true
                    try
                        nodeId=input('enter list of Node Id: ');
                    catch Error
                        disp(Error.message)
                        fprintf('Must specify a number or vector e.g. 1 or [1,2]\n')
                    end
                    %Check the type of input
                    if ~any(nodeId==0) && isvector(nodeId)
                        %Check that the input makes sense
                        flagBreak=true;
                        for j=1:length(nodeId) 
                           if and(~any(nodeId(j)==nId),enforceClusterOnlyMod)
                              flagBreak=false; 
                              fprintf('You specified a nodeId outside the cluster\n')
                           end
                        end
                        if flagBreak
                            break;
                        end
                    end
                end
                fid = fopen('nAddList.txt', 'a+');
                for j=1:length(nodeId)
                    fprintf(fid, '%d\n', nodeId(j));
                end
                fclose(fid);
            elseif isempty(option)
                fprintf('Moving on to the next grain\n');
                ind=ind+1;
                nextGrain=true;
            else
                %This should never happen.. but if it does catch it
                error('Error: unhandled grain')
            end
            if nextGrain
               break; 
            end
        end
        close(h);
    end
    fprintf('Finished all clusters in loop..\n');
    fprintf('Remember to run ClusterGrainsTwins for all edits to take place\n')
end
function G=removeEdge(G,rEdge,egroupId)
    rEdgeId=egroupId(rEdge);
    removeEdges=zeros(size(G.Edges.pairs,1),1,'logical');
    removeEdges(rEdgeId)=true;
    G=rmedge(G,G.Edges.pairs(removeEdges,1),...
        G.Edges.pairs(removeEdges,2));
end