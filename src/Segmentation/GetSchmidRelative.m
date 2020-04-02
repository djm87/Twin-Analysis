function [G_Family,edgeList] = GetSchmidRelative(G_Family,groupList,oriFamily,typeEdge,opt)
%The script returns the effective schmid factor on K1 plane in the eta1
%direction. Since this is called before we know the parent, we don't know 
%the sign of the twin variant rotation axis. For this reason variant 
%extraction is done in a seperate routine after the family tree is formed.

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
    sigmaPrinciple = sigma.eig;
    tauMax = (max(sigmaPrinciple) - min(sigmaPrinciple)) / 2;  

    %Allocate arrays for storage
    sigma13 = zeros(nedges,2);
    EffSF = zeros(nedges,2);

    %Loop over edges
    for i = 1:nedges
        type=typeEdgeReduced(i);
        %If the type is unknown then the effective schmid is zero and the
        %contribution in vote for schmid will be zero.
        if type <= nTwin && type > 0
            %Convert grains to s->c and apply symmetry operation
            gA = inv(grainsA(i).symmetrise); 
            gB = inv(grainsB(i).symmetrise);

            %Compute grain misorientation                
            mori = gA * grainsB(i).symmetrise;
            MM = mori(:,:) * twin{type}.RMT'; 

            %Get the smallest misorientation and its sym operation
            angs = angle(MM, 'noSymmetry') ./ degree;
            [~,id] = min(angs);
            [sym_ops1,sym_ops2] = ...
                ind2sub([sqrt(length(MM)) sqrt(length(MM))],id);

            %Get the stress in the twin frame i.e. s->c->t
            sigmaA = matrix(twin{type}.Rtw * (gA(sym_ops1) * sigma));
            sigmaB = matrix(twin{type}.Rtw * (gB(sym_ops2) * sigma));

            %Extract the stress on the k1 plane in eta1 direction
            sigma13(i,:) = [sigmaA(1,3),sigmaB(1,3)];

            %Compute the effective schmid factor
            EffSF(i,:) = [sigmaA(1,3) / (2*tauMax),sigmaB(1,3) / (2*tauMax)];
        end
    end %end loop over edges

    %Store arrays that were computed based on edgeList
    G_Family.Edges.sigma13(edgeList,:) = sigma13;
    G_Family.Edges.EffSF(edgeList,:) = EffSF;
end

