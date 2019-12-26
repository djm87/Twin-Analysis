function [] = UpdateClusters(G_clust,G_Initial,mergedGrains)
%UpdateClusters Recalculate the mergedGrains based on edges that have been
%removed or added. The UpdateClusters routine minimizes costly computation
%of family id and effSchmid

    eRemoveList=load('eRemoveList.txt'); %Remove the edge 
    eRemoveCleanup=load('eRemoveCleanup.txt'); %Edges determined by cleanFamilyTree
    eAddList=load('eAddList.txt'); %Adds edge and tries to give it a twin relationship, otherwise type twin unknown
    nRemoveList=load('nRemoveList.txt'); %Removes all edges associated with node
    nAddList=load('nAddList.txt'); %Tries to add all edges that have twin relationship
        
    
end

