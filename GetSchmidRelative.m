function G = GetSchmidRelative(G,twin,sigma)
%The script returns the effective schmid factor on K1 plane in the eta1
%direction. In addition, the symmetry operators are returned. Since this is
%called before we know the parent, we don't know the sign of the twin
%variant rotation axis. For this reason variant extraction is done in a
%seperate routine after the family tree is formed.

    % extract grains 
    grainIdA = G.Edges.pairs(:,1);
    grainIdB = G.Edges.pairs(:,2);
    grainsA = G.Nodes.meanOrientation(grainIdA);
    grainsB = G.Nodes.meanOrientation(grainIdB);
    
    %Compute the max shear with principle stresses
    %Remember that principle stresses are eigenvalues and the ordering of
    %the principle stress is from large (1) to small (3). Look at mohrs 
    %circle for info.
    sigmaPrinciple = sigma.eig;
    tauMax = (max(sigmaPrinciple) - min(sigmaPrinciple)) / 2;  
    
    %Get array sizes 
    ngrains = length(G.Nodes.Id);
    nedges = length(G.Edges.pairs);
    
    %Allocate arrays for storage
    sigma13 = zeros(nedges,2);
%     sym_ops = zeros(nedges,2);
    grainEffSF = zeros(ngrains, max(G.Edges.type));
  
    %Loop over edges
    for i = 1:nedges  
        
        type = G.Edges.type(i); 
        
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
        sigma13(i,1) = sigmaA(1,3);
        sigma13(i,2) = sigmaB(1,3);
                
        %Compute the effective schmid factor
        grainEffSF(grainIdA(i),type) = sigma13(i,1) / (2*tauMax);
        grainEffSF(grainIdB(i),type) = sigma13(i,2) / (2*tauMax);
        

    end %Loop over edges
    
    %Store arrays
    G.Edges.sigma13 = sigma13;
    G.Nodes.EffSF = grainEffSF;
    G.Edges.EffSFRelative = G.Edges.EffSF(:,1) - G.Edges.EffSF(:,2);
    
end

