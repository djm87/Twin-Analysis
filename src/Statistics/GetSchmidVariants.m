function [schmidVariantRank,k1NormalAngle]=GetSchmidVariants(G_Family_sub,nParentFamilyId,nType,opt)
%The script returns the effective schmid factor on K1 plane in the eta1
%direction for all variants along with the active variant

    %Find twinned grains 
    isTwin=nParentFamilyId>0;
    grainsC=G_Family_sub.Nodes.meanOri(isTwin);
    grainsP=G_Family_sub.Nodes.meanOri(nParentFamilyId(isTwin));   
    typeC=nType(isTwin);
    
    %Compute the max shear with principle stresses
    %Remember that principle stresses are eigenvalues and the ordering of
    %the principle stress is from large (1) to small (3). Look at mohrs 
    %circle for info.
%     sigmaPrinciple = opt.sigma.eig;
%     tauMax = (max(sigmaPrinciple) - min(sigmaPrinciple)) / 2;  
    
    %Get table array sizes 
    ngrains = length(grainsP);
    
    %Get the maximum number of variants over all twin modes
    maxVariants = 0;
    for i=1:opt.nTwin
        tmp = length(opt.twin{i}.axisVariants);
        if tmp > maxVariants
           maxVariants = tmp; 
        end
    end
    
    %Allocate arrays for tmp storage
    k1NormalAngle_tmp = zeros(ngrains, 1);
    schmidVariantRank_tmp = zeros(ngrains, 1);
    
    %Loop over edges
    for i = 1:ngrains  
        
        %For each possible variant determine the ideal child grain and the
        %corresponding schmid factor
        %Get the rotation axes
        vari = grainsP(i)*opt.twin{typeC(i)}.axisVariants;

        %set a rotation around that axis
        rot = rotation('axis',vari,'angle',opt.twin{typeC(i)}.angle*degree);

        % rotate the parent around the twin axis,
        oriV = rot*grainsP(i);        
        mis = angle(oriV, grainsC(i)) / degree; %all the varients of the twin mode

        %determines active twin variant
        [~,activeVariant] = min(abs(mis)); 

        %determine the schmid on each twin variant. The schmid in the twin
        %is the negative of the parent. We want the parent stats so multipy
        %by -
        schmidVariants = -opt.twin{typeC(i)}.sS.SchmidFactor( oriV \ opt.sigma);

        %rank the twin variant
        [~,schmidVariantRank] = sort(schmidVariants,'descend');

        schmidVariantRank_tmp(i) = find(schmidVariantRank==activeVariant);
        if length(opt.twin{typeC(i)}.k1)==1
            k1NormalAngle_tmp(i)=cos(angle(grainsC(i) \ vector3d(opt.twin{typeC(i)}.k1),zvector)); %Make sure rot is applied correctly
        end

    end %end main loop
    
    k1NormalAngle_tmp(k1NormalAngle_tmp<0)=k1NormalAngle_tmp(k1NormalAngle_tmp<0)+pi;
    k1NormalAngle=zeros(length(nParentFamilyId),1);
    k1NormalAngle(isTwin)=k1NormalAngle_tmp;
    
    schmidVariantRank=zeros(length(nParentFamilyId),1);
    schmidVariantRank(isTwin) = schmidVariantRank_tmp;
    
end

