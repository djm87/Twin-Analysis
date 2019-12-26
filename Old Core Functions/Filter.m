function [ePARISChange]=RemovalPARISImprovement(G,grains,singleFragRelationship,twinModes)
%Filter currently just performs PARIS based edge removal
    [mergedGrains,parentId,combinedTwinBoundary,combine] = MergeByEdge(G.Edges.pairs,G.Edges.combineCleaned,grains);
    paris_orig=mergedGrains.paris;
    %     figure;
    %     plot(grains,grains.meanOrientation);hold on;
    %     plot(mergedGrains.boundary,'linecolor','k','linewidth',2,'linestyle','-',...
    %         'displayName','merged grains'); hold off;
    toCombine=G.Edges.combineCleaned;
    ePARISChange=zeros(length(G.Edges.combineCleaned),1);
    %loop over groups and identify edges to try removing
    isInside = checkInside(grains, grains);
    [GrainIdInclusion,~] = find(isInside);
    
    for i=1:max(G.Edges.Group)
        %load edge and node properties for clustered fragments
        egroupId = find((i==G.Edges.Group)==true); %converts logical arrays to indices
        ngroupId = find((i==G.Nodes.Group)==true);
        nFamily = G.Nodes.FamilyID(ngroupId);
        eType = G.Edges.type(egroupId);
        ePairs = G.Edges.pairs(egroupId,:);
        eFamily = G.Edges.FamilyID(egroupId,:);
        
        %We won't consider removing a family with multiple fragments
        [familyCnt,family]=hist(nFamily,unique(nFamily));
        families2Remove=family(familyCnt==1);
        nFamilies2Remove=length(families2Remove);
        
        %determine which edges have the family and try removing (based on
        %filter options)
        for j=1:nFamilies2Remove
            e2Remove=find(sum(eFamily==families2Remove(j),2));
            if length(e2Remove)==1 && singleFragRelationship
                n1=ePairs(e2Remove,1);
                n2=ePairs(e2Remove,2);
                
                %We should ignore internal grains and edge types that
                %are excluded
                if ~any(n1==GrainIdInclusion | n2==GrainIdInclusion) && any(eType(e2Remove)==twinModes)
                    if size(toCombine,2)<j
                        toCombine(:,j)=G.Edges.combineCleaned;
                    end
                    toCombine(egroupId(e2Remove),j)=false;
                    
                end
            elseif ~tFilter.singleFragRelationship
                %Find if all edges can be removed
                canRemoveTypes=true;
                for k=1:length(e2Remove)
                    if ~any(eType(e2Remove(k))==tFilter.twinModes)
                        canRemoveTypes=false;
                    end
                end
                for k=1:length(e2Remove)
                    n1=ePairs(e2Remove(k),1);
                    n2=ePairs(e2Remove(k),2);
                    
                    %We should ignore internal grains and edge types that
                    %are excluded
                    if ~any(n1==GrainIdInclusion | n2==GrainIdInclusion) && canRemoveTypes
                        if size(toCombine,2)<j
                            toCombine(:,j)=G.Edges.combineCleaned;
                        end
                        toCombine(egroupId(e2Remove(k)),j)=false;
                    end
                end
            end
        end
    end
    
    %for each edge(s) removed determine if it improved anything
    for cnt=1:size(toCombine,2)
        toCombine_local=toCombine(:,cnt);
        [mergedGrains2,parentId2,~,~] = MergeByEdge(G.Edges.pairs,toCombine_local,grains);
        paris_split=mergedGrains2.paris;
        for i=1:max(G.Edges.Group)
            %load edge and node properties for clustered fragments
            egroupId = find((i==G.Edges.Group)==true); %converts logical arrays to indices
            ngroupId = find((i==G.Nodes.Group)==true);
            nId = G.Nodes.Id(ngroupId);
            
            %first find if we tried removing any edges
            if any(~toCombine_local(egroupId))
                %the nodes we are interested in
                nId_orig=unique(mergedGrains(parentId(nId)).id);
                nId_split=mergedGrains2(parentId2(nId)).id;
                nId_splitu=unique(nId_split);
                
                %If we break a frag from a cluster only look at cluster for
                %the paris change eval
                if length(nId_splitu)<length(nId_split)
                    [nIdCnt,nId_split]=hist(nId_split,nId_splitu);
                    nId_split=nId_split(nIdCnt>1);
                end
                area=mergedGrains2(nId_split).area;
                mean_paris_split=dot(area,paris_split(nId_split))/sum(area);
                ePARISChange(~toCombine_local(egroupId))=(paris_orig(nId_orig)-mean_paris_split)/paris_orig(nId_orig)
            end
        end
    end
end

