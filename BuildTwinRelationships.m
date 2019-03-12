function [G_Complete] = BuildTwinRelationships(G_clust,grains,w,twin,sigma,seg_angle_grouped,Mistol,meanMistolRelaxed)
% This function computes votes, cleans, tree, and builds the twin map from 
%which stats can be computed  Detailed explanation goes here
    %% Compute Schmid info for twin/parents in clustered grains
    %This computes Schmid factor for twin/parent identification

    %Should be updated with the mtex method as in GetSChmidVariants
    G_clust = GetSchmidRelative(G_clust,twin,sigma);

    tmp=G_clust.Nodes.EffSF(:,1);
    notZero=tmp<0;
      figure; 
        plot(grains(G_clust.Nodes.Id(notZero)),...
            tmp(notZero),'Micronbar','off');
        hold on 
        p=plot(G_clust,'XData',G_clust.Nodes.centroids(:,1),...
            'YData',G_clust.Nodes.centroids(:,2));
        hold off
    drawnow
    %% Identify Families in grain clusters 
    %Remove nodes that aren't connected by edges, alternatively use condition
    %on minimum occurance from conncomp
    doplot=true;
    dolabel=true;
    G_clust = AssignFamilyIDs(G_clust,grains,seg_angle_grouped,doplot,dolabel);
    drawnow

    %% Implement voting scheme for the respective families
    % For each family compute the vote between connected families 
    %Edges currently connect clustered grains. We can compute the vote for each
    %edge and then average for the family. This ensures that only neighbors are
    %ever voted on. So that votes are really between familys average values of
    %schmid are computed and stored for each fragment in a family. The area and
    %relative boundary are taken in the addative sense rather than average. 

    %Edges.Gb contains boundaries between pairs 
    %Nodes.Gb contains boundaries for each fragment 
    %To compute the relative boundary between families we need to sum boundary
    %of a particular type and since each pair is its own type this should be
    %relatively simple. 
    G_Complete_unclean = FamilyVotes(G_clust,w);

    %% Cleanup the family tree
    runCleanup=true
    maxIter=5;
    CleanupIter=0;
    %I think this shouldn't take more than 2 itterations. There can be some
    %cases that aren't solved until grains are regroupd (i.e you have two
    %seperate grains that need still need to be split and a circular relationship exists
    %because the same families are in both grains. It's a thing.. go figure!)
    while runCleanup && CleanupIter<maxIter
        [G_Complete_unclean,runCleanup] = CleanFamilyTree(G_Complete_unclean,grains);

        if runCleanup
            G_Complete_unclean=ClusterGrainsTwins(G_Complete_unclean,grains,Mistol,meanMistolRelaxed,...
            twin,[],false,false);
        end

        CleanupIter=CleanupIter+1;
    end
    if CleanupIter == maxIter
        disp('More CleanupIter than expected, likely need to solve something manually or improve the code')
    end
    %% Reconstruct grains based on cleanup (some new grains may have been created
    edgeList2Remove=[]; % number of edge from label 

    doplot=true;
    dolabel=false;
    G_Complete = ClusterGrainsTwins(G_Complete_unclean,grains,Mistol,meanMistolRelaxed,...
        twin,edgeList2Remove,doplot,dolabel)
    drawnow
    %% Redo Families in grain clusters 
    %Remove nodes that aren't connected by edges, alternatively use condition
    %on minimum occurance from conncomp
    doplot=true;
    dolabel=false;
    G_Complete = AssignFamilyIDs(G_Complete,grains,seg_angle_grouped,doplot,dolabel);
    drawnow
    %% Create the family tree 
    %uses the recursive relationship as described in Pradalier et al. 2018
    [G_Complete] = CreateFamilyTree(G_Complete,grains);
end

