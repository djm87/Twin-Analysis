function [G_Complete] = CreateFamilyTree(G_Complete,grains)
%CreateFamilyTree takes a clean set of family/edge relationships
%and assigns a twin type and generation to each fragment.
%
%Note: Garbage in, garbage out. There are no checks in this routine that 
%makes sure things are clean and make sense.
    
    %Intialize new output arrays
    G_Complete.Nodes.Generation = zeros(length(G_Complete.Nodes.Id),1,'int8');
    G_Complete.Nodes.Type = zeros(length(G_Complete.Nodes.Id),1,'int8');
    G_Complete.Nodes.TypeColored = zeros(length(G_Complete.Nodes.Id),3,'int8');
    colors=hsv(double(max(G_Complete.Edges.type))+1);
    %Loop over the groups 
    for i=1:max(G_Complete.Edges.Group) 
        egroupId = find((i==G_Complete.Edges.Group)==true); %converts logical arrays to indices
        ngroupId = find((i==G_Complete.Nodes.Group)==true);
        nId = G_Complete.Nodes.Id(ngroupId);
        eType = G_Complete.Edges.type(egroupId);
        ePairs = G_Complete.Edges.pairs(egroupId,:);
        eParent = G_Complete.Edges.Parent(egroupId,:);
        nGeneration = -ones(length(nId),1,'int8');
        nType = -ones(length(nId),1,'int8');
        
        %Get the edge connectivity matrix
        EdgeMatrix=zeros(length(nId),length(nId),'uint8');
        for k = 1:size(ePairs,1)
            p=find(ePairs(k,eParent(k,:))==nId);
            c=find(ePairs(k,~eParent(k,:))==nId);
            EdgeMatrix(p,c)=eType(k); %Type and relationship (by position)
        end
        
        %Call labeling function
        
        [nGeneration,nType] = fillGenerations(EdgeMatrix,...
                                            nGeneration,nType,[],0,8);
        G_Complete.Nodes.Type(ngroupId)=nType;
        G_Complete.Nodes.Generation(ngroupId)=nGeneration;
        G_Complete.Nodes.TypeColored(ngroupId,:)=colors(nType+1,:);

    end
end
function [nGeneration,nType] = fillGenerations(EdgeMatrix,...
                                        nGeneration,nType,genId,gen,maxGen)
    %Simple recursive function to label type and generation 
    if gen==0
        genId = find(sum(EdgeMatrix,1)==0);
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
            genId1 = find(EdgeMatrix(genId(i),:)~=0);
            if ~isempty(genId1)
                nGeneration(genId1)=gen;
                nType(genId1)=EdgeMatrix(genId(i),genId1);
                if any(nGeneration==-1)
                    [nGeneration,nType] = fillGenerations(EdgeMatrix,...
                        nGeneration,nType,genId1,gen+1,maxGen);
                end
            end
        end
        if any(nGeneration==-1) && gen==1
            nGeneration(:)=-1;
            nType(:)=0;
            disp('Warning: Failure to relate all fragments')
        end
    end
    
end