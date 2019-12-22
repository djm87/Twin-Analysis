function [G,err_group] = CreateFamilyTree(G,grains,mergedGrains,twin)
%CreateFamilyTree takes a clean set of family/edge relationships
%and assigns a twin type and generation to each fragment.
%
%Note: Garbage in, garbage out. There are no checks in this routine that 
%makes sure things are clean and make sense. Clean family tree should be
%used to remove redundant relationships and make sure things are good.
    
    %Intialize new output arrays
    G.Nodes.Generation = zeros(length(G.Nodes.Id),1,'int8');
    G.Nodes.Type = zeros(length(G.Nodes.Id),1,'int8');
    G.Nodes.TypeColored = zeros(length(G.Nodes.Id),3,'int8');
    G.Nodes.noParent = zeros(size(G.Nodes.Id,1),1,'logical');
    G.Edges.noParent = zeros(size(G.Edges.pairs,1),1,'logical');
    colors=hsv(double(max(G.Edges.type))+1);
    
    %Loop over the groups 
    openType=length(twin);
    groups=unique(G.Edges.Group);
    err_group=zeros(max(groups),1,'logical');
    for i=1:length(groups)
        group=groups(i);
        %load edge and node properties for clustered fragments
        egroupId = find((group==G.Edges.Group)==true); %converts logical arrays to indices
        ngroupId = find((group==G.Nodes.Group)==true);
        typeKnown=~G.Nodes.typeUnknown(ngroupId);
        nId = G.Nodes.Id(ngroupId(typeKnown));
        nFamily = G.Nodes.FamilyID(ngroupId(typeKnown));
        nGeneration = -ones(length(nId),1,'int8');
        nType = -ones(length(nId),1,'int8');        
        eType = G.Edges.type(egroupId);
        etypeKnown=eType~=openType;
        eType=eType(etypeKnown);
        egroupId=egroupId(etypeKnown);
        eParent = G.Edges.Parent(egroupId,:);
        ePairs = G.Edges.pairs(egroupId,:);        
        eFamily = G.Edges.FamilyID(egroupId,:);


        %Check to see if some of the parent relationships were not figured
        %out
%         noParent=all(eParent==0,2);
%         if any(noParent)
%             G.Edges.noParent(egroupId(noParent))=true;
%             nId_noParent=unique(ePairs([noParent,noParent]));
%             egroupId=egroupId(~noParent);
%             eType = G.Edges.type(egroupId);
%             eParent = G.Edges.Parent(egroupId,:);
%             ePairs = G.Edges.pairs(egroupId,:);
%             nId=unique(ePairs);
%             for j=1:length(nId_noParent)
%                if all(nId_noParent(j)~=nId)
%                    G.Nodes.noParent(nId_noParent(j))=true;
%                end
%             end
%         end
        
%         figure; plot(grains(nId),G.Nodes.FamilyID(nId))
        if ~isempty(eType)       
            %Get the edge connectivity matrix
            EdgeMatrix=zeros(length(nId),length(nId),'uint8');
            for k = 1:size(ePairs,1)
                p=find(ePairs(k,eParent(k,:))==nId);
                c=find(ePairs(k,~eParent(k,:))==nId);
                EdgeMatrix(p,c)=eType(k); %Type and relationship (by position)
            end
            if(size(EdgeMatrix,1)~=size(EdgeMatrix,2))
                error('EdgeMatrix not square')
            end
            %Call labeling function
            [nGeneration,nType,err] = fillGenerations(EdgeMatrix,...
                                                nGeneration,nType,[],0,8,false);
            
            %catch any issues and flag the group
            if err 
                err_group(group)=true;
                FamilyMatrix = buildFamilyMatrix(nFamily,eFamily,eParent)
                value=G.Nodes.meanOrientation;
                plotNeighbors=true;
                enforceClusterOnlyMod=false;
                [G,~,~] = ClusterEditor(group,G,grains,mergedGrains,value,0,plotNeighbors,enforceClusterOnlyMod); 
            end
            
            %Set the twin based properties
            G.Nodes.Type(ngroupId(typeKnown))=nType;
            G.Nodes.Generation(ngroupId(typeKnown))=nGeneration;
            G.Nodes.TypeColored(ngroupId(typeKnown),:)=colors(nType+1,:);

            %Set the unknown grain properties 
            if any(~typeKnown)
                G.Nodes.Type(ngroupId(~typeKnown))=openType;
                G.Nodes.Generation(ngroupId(~typeKnown))=-1;
                G.Nodes.TypeColored(ngroupId(~typeKnown),:)=colors(G.Nodes.Type(ngroupId(~typeKnown))+1,:);
            end
        end
    end
end
function [nGeneration,nType,err] = fillGenerations(EdgeMatrix,...
                                        nGeneration,nType,genId,gen,maxGen,err)
    %Simple recursive function to label type and generation 
    if gen==0
        genId = find(all(EdgeMatrix==0,1)); %Columns are children
        if ~isempty(genId)
            nType(genId)=0;
            nGeneration(genId)=0;
            gen=1;
        else
            nType(:)=0;
            nGeneration(:)=-1;
            disp('Warning: No parent generation is obvious')
        end
    end
    if gen==maxGen
        %Reset nGeneration so we know to debug.
        nGeneration(:)=-1;
        nType(:)=0;
    else
        for i=1:length(genId)
%             i
%             EdgeMatrix
%             genId
            genId1 = find(EdgeMatrix(genId(i),:)~=0);
            if ~isempty(genId1)
                nGeneration(genId1)=gen;
                nType(genId1)=EdgeMatrix(genId(i),genId1);
                if any(nGeneration==-1)
                    [nGeneration,nType,err] = fillGenerations(EdgeMatrix,...
                        nGeneration,nType,genId1,gen+1,maxGen,err);
                end
            end
        end
        if any(nGeneration==-1) && gen==1
            nGeneration(:)=-1;
            nType(:)=0;
            disp('Warning: Failure to relate all fragments')
            err=true;
        end
    end
    
end

function FamilyMatrix = buildFamilyMatrix(nFamily,eFamily, Parent)
    FamilyMatrix=zeros(max(nFamily),'logical');
    %Make Family Tree 
    for j = 1:size(eFamily,1)
        p=eFamily(j,Parent(j,:));
        c=eFamily(j,~Parent(j,:)); 
        FamilyMatrix(p,c)=true;
    end
end