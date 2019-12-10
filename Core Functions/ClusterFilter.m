function [G,runCleanupAgain,ind] = ClusterEditor(groups,G,grains,value,ind)
%ClusterFilter loops through clusters in prefiltered groups, allowing a 
%user to add/remove nodes edges and relationships
%Inputs:  groups - list of merged grain Ids 
%         G - graph
%         grains - 2dgrain object
%         value - matches length grains and contains a scalar value or orientation
%Outputs: Most outputs go to .txt files used in ClusterGrains.
%         The Graph object is updated if only if an edge is deleted.
%         runCleanupAgain is used in CleanFamilyTree to recluster grains
    rEdge=[];
    runCleanupAgain=[];
    for i=1:length(groups) 
        group=groups(i);
        egroupId = find((group==G.Edges.Group)==true); %converts logical arrays to indices
        ngroupId = find((group==G.Nodes.Group)==true);
        nId = G.Nodes.Id(ngroupId);
        nFamily = G.Nodes.FamilyID(ngroupId);
        eGlobalId = G.Edges.GlobalID(egroupId);  
        
        %Plot grain cluster
        h=figure; plot(grains(nId),value(nId))
%         text(grains(nId),int2str(nId))
        hold on
        toremove=ones(length(G.Nodes.Id),1,'logical');
        toremove(unique(G.Edges.pairs(egroupId,:)))=false;
        G_Removed=rmnode(G,find(toremove));
        p=plot(G_Removed,'XData',G_Removed.Nodes.centroids(:,1),...
            'YData',G_Removed.Nodes.centroids(:,2),'displayName','graph');
        pairs1=G_Removed.Edges.pairs(:,1);
        pairs2=G_Removed.Edges.pairs(:,2);
        for j=1:length(G_Removed.Nodes.Id)
            pairs1(pairs1==G_Removed.Nodes.Id(j))=j;
            pairs2(pairs2==G_Removed.Nodes.Id(j))=j;
        end
        labeledge(p,pairs1,...
            pairs2,G_Removed.Edges.GlobalID); 
        hold off   

        %give node and edge info
        fprintf('===================================\n')
        fprintf('Group %d\n',i)
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
        fprintf('1 to specify node as twin and parent of none\n')
        fprintf('2 to specify node as parent of none\n')
        fprintf('3 to set an edge relationship\n')
        fprintf('4 to remove an edge\n')
        
        fprintf('Any other input will abort\n')
        option=input('enter number: ');
        if option==1
            nodeId=input('enter list of Node Id: ');
            fid = fopen('notParent.txt', 'a+');
            for j=1:length(nodeId)
                fprintf(fid, '%d\n', nodeId(j));
            end
            fclose(fid);
        elseif option == 2
            nodeId=input('enter list of Node Id: ');
            fid = fopen('notTwin.txt', 'a+');
            for j=1:length(nodeId)
                fprintf(fid, '%d\n', nodeId(j));
            end
            fclose(fid);
        elseif option == 3  
            fprintf('To specify a relation enter edge id, and logical pair e.g. 124 1 0 where 1 is the parent\n')
            promt=='y';
            while promt=='y'
                eRelation=input('specify edge relation: ');
                fid = fopen('eRelation.txt', 'a+');
                fprintf(fid, '%d %d %d\n', eRelation);
                fclose(fid);
                prompt=input('Would you like to add another relation? (y or n): ','s');
            end
        elseif option == 4                    
            edgeId=input('enter list of edge Id: ');
            fid = fopen('eRemoveList.txt', 'a+');
            for j=1:length(edgeId)
                fprintf(fid, '%d\n', edgeId(j));
            end
            fclose(fid);
            rEdge=find(edgeId==eGlobalId);
            G=removeEdge(G,rEdge,egroupId);
            ind=ind+1;
            runCleanupAgain=true;            
        elseif isempty(option)
            %Move on to the next grain
            ind=ind+1;
        else
            error('Error: unhandled grain')
        end
        close(h);
    end
end