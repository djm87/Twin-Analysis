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
    sigmaPrinciple = sigma.eig;
    tauMax = (max(sigmaPrinciple) - min(sigmaPrinciple)) / 2;  

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
%             SFVA = twin{type}.sS.SchmidFactor(grainsA(i) \ sigma)
            aVA = activeSchmidVariant(grainsA(i),grainsB(i),twin{type}.axisVariants,twin{type}.angle);
            
%             SFVB = twin{type}.sS.SchmidFactor(grainsB(i) \ sigma)
            aVB = activeSchmidVariant(grainsB(i),grainsA(i),twin{type}.axisVariants,twin{type}.angle);
%             EffSF(i,:) = [SFVA(aVA),SFVB(aVB)];
            
            sigmaA=matrix(twin{type}.Rtw * inv(grainsA(i)) * sigma);
            sigmaB=matrix(twin{type}.Rtw * inv(grainsB(i)) * sigma);
            reshape(sigmaA(1,3,:),6,1)

            %Extract the stress on the k1 plane in eta1 direction
            sigma13(i,:) = [sigmaA(1,3,aVB),sigmaB(1,3,aVA)];

            %Compute the effective schmid factor
            EffSF(i,:) = [sigmaA(1,3,aVA) / (2*tauMax),sigmaB(1,3,AVB) / (2*tauMax)];
            
            %Active variant 
            SFAV(i,:) = [aVA,aVB];
            
            %Variant Rank 
            [~,SFVRA] = sort(sigmaA(1,3,:),'descend');
            [~,SFVRB] = sort(sigmaB(1,3,:),'descend');
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

function activeVariant=activeSchmidVariant(gA,gB,axisVariants,theta)
    %The active variant index is with respect to a fixed set of twin variants and corresponds to the
    %twin slip systems. Therefore the effective schmid is also available.
    
    %Rotate the variants to the sample frame of gA
    vari = gA*axisVariants;

    %set a rotation around those axes 
    rot = rotation('axis',vari,'angle',theta*degree);

    % rotate the parent around the twin axis,
    oriV = rot*gA;
    
    %compute the misorientation between twinned gA orientations and the parent
    mis = angle(oriV, gB) / degree; %all the varients of the twin mode

    %the active twin variant for a grain determined to be twinning is the minimum
    [~,activeVariant] = min(abs(mis)); 
end