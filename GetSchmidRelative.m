function G = GetSchmidRelative(G,twin,sigma)
%The script returns the effective schmid factor on K1 plane in the eta1
%direction. In addition, the symmetry operators are returned so the exact
%twin variant is known.

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
    
    %Allocate arrays for storage
    ngrains = length(G.Nodes.Id);
    nedges = length(G.Edges.pairs);
    sigma13 = zeros(nedges,2);
    sym_ops = zeros(nedges,2);
    grainEffSF = zeros(ngrains, max(G.Edges.type));
    
    %Loop over twin types
    for type = 1:length(twin)
        %Set the symmetry operator for the twin
        CS = twin{type}.CS;
    
        %Vectorize what we can for grain misorientation
%         Lg2=CS(:) * grainsB;
%         g1g2tLt = Lg2;
%         for i = 1:length(CS)
%             g1g2tLt(i,:) = grainsA .* transpose(Lg2(i,:)');
%         end
        
        %Loop over edges
        for i = 1:nedges  
            if G.Edges.type(i) == type
                %Compute grain misorientation
                Lg2=CS(:)*grainsB(i);
                g1g2tLt=(grainsA(i)*Lg2');

                M = CS(:) * g1g2tLt;

                %Find the difference between twin and grain misorientation
                MM = M * twin{type}.RMT;
                angs = angle(MM, 'noSymmetry')' / degree;

                %Extract the closest match along with the symmetry operations
                [~,id] = min(angs);
                [sym_ops(i,1),sym_ops(i,2)] = ...
                    ind2sub([length(CS) length(CS)],id);

                %Find the resolved shear stress on the K1 plane in the eta1
                %direction
                aA=twin{type}.Rtw.matrix * CS(sym_ops(i,1)).matrix * grainsA(i).matrix;
                aB=twin{type}.Rtw.matrix * CS(sym_ops(i,2)).matrix * grainsB(i).matrix;

                sigma13(i,1) = dot(aA(1,:), aA(3,:) * sigma.matrix');
                sigma13(i,2) = dot(aB(1,:), aB(3,:) * sigma.matrix');

                grainEffSF(grainIdA(i),type) = sigma13(i,1) / (2*tauMax);
                grainEffSF(grainIdB(i),type) = sigma13(i,2) / (2*tauMax);
            end
        end %Loop over edges
    end %Loop over twin types
    
    %Store arrays
    G.Nodes.EffSF = grainEffSF;
    G.Edges.sigma13 = sigma13;
    G.Edges.EffSF = sigma13 ./ (2 * tauMax);
    G.Edges.EffSFRelative = G.Edges.EffSF(:,1) - G.Edges.EffSF(:,2);
    
end

