function [G_Family,edgeList] = GetSchmidRelative(G_Family,groupList,oriFamily,typeEdge,opt)
%The script returns the effective schmid factor on K1 plane in the eta1
%direction.

    %Find the subset of edges to do computation for
    edgeList=intersect_wrepeat(groupList,G_Family.Edges.Group);
    
    %Exit routine if there is nothing to do
    if isempty(edgeList)
       return 
    end

    %Get array sizes 
    nedges = length(edgeList);

    % extract mean orientation for each family edge relationship in group
    grainsA = oriFamily(edgeList,1);
    grainsB = oriFamily(edgeList,2);

    % extract options 
    sigma=opt.sigma; 
    twin=opt.twin;
    nTwin=opt.nTwin;
    typeEdgeReduced=typeEdge(edgeList);
    %Compute the max shear with principle stresses
    %Remember that principle stresses are eigenvalues and the ordering of
    %the principle stress is from large (1) to small (3). Look at mohrs 
    %circle for info.
%     sigma=stressTensor([1 1 0; 0 -1 0; 0 0 0])
    %sigmaPrinciple = sigma.eig;
    %tauMax = (max(sigmaPrinciple) - min(sigmaPrinciple)) / 2;  

    %Allocate arrays for storage
    sigma13 = zeros(nedges,2);
    EffSF = zeros(nedges,2);
    SFAV = zeros(nedges,2);
    SFAVR = zeros(nedges,2);
    
    %Loop over edges
    for i = 1:nedges
        type=typeEdgeReduced(i);
        %If the type is unknown then the effective schmid is zero and the
        %contribution in vote for schmid will be zero.
        if type <= nTwin && type > 0
            %Convert grains to s->c and apply symmetry operation
            %sSTwin= grainsA(i) * twin{type}.sS
            %SFVA = sSTwin.SchmidFactor(sigma)
            SFVA = twin{type}.sS.SchmidFactor(inv(grainsA(i)) * sigma);
            aVA = activeSchmidVariant(grainsA(i),grainsB(i),twin{type}.RMTall);

            SFVB = twin{type}.sS.SchmidFactor(inv(grainsB(i)) * sigma);
            aVB = activeSchmidVariant(grainsB(i),grainsA(i),twin{type}.RMTall);
            EffSF(i,:) = [SFVA(aVA),SFVB(aVB)];
            
            
            %sigmaA=matrix(twin{type}.Rtw * inv(grainsA(i)) * sigma);
            %sigmaB=matrix(twin{type}.Rtw * inv(grainsB(i)) * sigma);
            %reshape(sigmaA(1,3,:),12,1)
            %reshape(sigmaB(1,3,:),12,1)
            
            %sS = slipSystem(twin{i}.CS.cAxis,twin{i}.CS.aAxis)
            %SFVA = sS.SchmidFactor(twin{type}.Rtw * inv(grainsA(i)) * sigma)

            
            %Extract the stress on the k1 plane in eta1 direction
            %sigma13(i,:) = [sigmaA(1,3,aVA),sigmaB(1,3,aVB)];

            %Compute the effective schmid factor
            %EffSF(i,:) = [sigmaA(1,3,aVA) / (2*tauMax),sigmaB(1,3,aVB) / (2*tauMax)];
            
            %Active variant 
            SFAV(i,:) = [aVA,aVB];
            
            %Variant Rank 
            [~,SFVRA] = sort(SFVA,'descend');
            [~,SFVRB] = sort(SFVB,'descend');
            SFAVR(i,1) = find(SFVRA==aVA);            
            SFAVR(i,2) = find(SFVRB==aVB);
            
        end
    end %end loop over edges

    %Store arrays that were computed based on edgeList
    G_Family.Edges.sigma13(edgeList,:) = sigma13;
    G_Family.Edges.EffSF(edgeList,:) = EffSF;
    G_Family.Edges.SFAV(edgeList,:) = SFAV;    
    G_Family.Edges.SFAVR(edgeList,:) = SFAVR;

end

function childId=activeSchmidVariant(parentOri,childOri,p2t)
    %The active variant index is with respect to a fixed set of twin variants and corresponds to the
    %twin slip systems. Therefore the effective schmid is also available.
    
    % all child variants
    childVariants  = parentOri * inv(p2t);

    if size(childVariants,1) == 1
      childVariants = repmat(childVariants,length(childOri),1);
    end

    % compute distance to all possible variants
    d = dot(childVariants,repmat(childOri(:),1,size(childVariants,2)));

    % take the best fit
    [~,childId] = max(d,[],2);
end