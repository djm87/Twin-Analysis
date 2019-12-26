function G = GetSchmidRelative(G,twin,sigma)
%The script returns the effective schmid factor on K1 plane in the eta1
%direction for all variants along with the active variant

    %Get array sizes 
    ngrains = length(G.Nodes.Id);
    nedges = length(G.Edges.pairs);   
    openType = length(twin);
    
    % extract grain id and grains 
    grainIdP=ones(nedges,1);
    grainIdC=ones(nedges,1);
    for i=1:nedges
        if G.Edges.type(i) ~= openType
            grainIdP(i) = G.Edges.pairs(i,G.Edges.Parent(i,:));
            grainIdC(i) = G.Edges.pairs(i,~G.Edges.Parent(i,:));
        end
    end
    grainsP = G.Nodes.meanOrientation(grainIdP);
    grainsC = G.Nodes.meanOrientation(grainIdC);
    
    
    %Compute the max shear with principle stresses
    %Remember that principle stresses are eigenvalues and the ordering of
    %the principle stress is from large (1) to small (3). Look at mohrs 
    %circle for info.
    sigmaPrinciple = sigma.eig;
    tauMax = (max(sigmaPrinciple) - min(sigmaPrinciple)) / 2;  
    
    %Get table array sizes 
    ngrains = length(G.Nodes.Id);
    nedges = length(G.Edges.pairs);
    
    %Get the maximum number of variants over all twin modes
    maxVariants = 0;
    for i=1:length(twin)-1
        tmp = length(twin{i}.axisVariants);
        if tmp > maxVariants
           maxVariants = tmp; 
        end
    end
    
   
    %Allocate arrays for tmp storage
    schmid = cell(nedges, 1);
    schmidRank = cell(nedges, 1); 
    schmidActive = zeros(nedges, 1);
    schmidActiveRank = zeros(nedges, 1);
    schmidActiveN = zeros(nedges, 1);
    k1NormalAngle_tmp = zeros(nedges, 1);
    
    typeEdge = G.Edges.type; 
    %Loop over edges
    parfor i = 1:nedges  
        
        type = typeEdge(i); 
        if type ~= length(twin)
            %For each possible variant determine the ideal child grain and the
            %corresponding schmid factor
            %Get the rotation axes
            vari = grainsP(i)*twin{type}.axisVariants;

            %set a rotation around that axis
            rot = rotation('axis',vari,'angle',twin{type}.angle*degree);

            % rotate the parent around the twin axis,
            oriV = rot*grainsP(i);        
            mis = angle(oriV, grainsC(i)) / degree; %all the varients of the twin mode

            %determines active twin variant
            [~,activeVariant] = min(abs(mis)); 

            %determine the schmid on each twin variant. The schmid in the twin
            %is the negative of the parent. We want the parent stats so multipy
            %by -
            schmidVariants = -twin{type}.sS.SchmidFactor( oriV \ sigma);

            %rank the twin variant
            [~,schmidVariantRank] = sort(schmidVariants,'descend');

            schmid{i} = schmidVariants;
            schmidRank{i} = schmidVariantRank;
            schmidActiveN(i) = activeVariant;
            schmidActiveRank(i) = find(schmidVariantRank==activeVariant);
            schmidActive(i) = schmidVariants(activeVariant);
            if length(twin{type}.k1)==1
                k1NormalAngle_tmp(i)=cos(angle(oriV(activeVariant) \ vector3d(twin{type}.k1),zvector)); %Make sure rot is applied correctly
            else
                k1NormalAngle_tmp(i)=0
            end
        end
    end %Loop over edges
    k1NormalAngle = zeros(ngrains,1);
    k1NormalAngle(grainIdC)=k1NormalAngle_tmp;
    
    %Store arrays in table
    G.Edges.schmid = schmid;
    G.Edges.schmidRank = schmidRank;
    G.Edges.schmidActive = schmidActive;
    G.Edges.schmidActiveRank = schmidActiveRank;
    G.Edges.schmidActiveN = schmidActiveN;
    G.Nodes.k1NormalAngle = k1NormalAngle;
end

