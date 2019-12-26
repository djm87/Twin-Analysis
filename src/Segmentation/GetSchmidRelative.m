function [G_clust,G] = GetSchmidRelative(G_clust,G,grains,mergedGrains,opt)
%The script returns the effective schmid factor on K1 plane in the eta1
%direction. In addition, the symmetry operators are returned. Since this is
%called before we know the parent, we don't know the sign of the twin
%variant rotation axis. For this reason variant extraction is done in a
%seperate routine after the family tree is formed.
    %Only start calculations if there are some to be done
    doCalc=~G_clust.Edges.SFCalcd;
    if ~any(doCalc)
        return
    end
    
    % extract grains 
    grainIdA = G_clust.Edges.pairs(:,1);
    grainIdB = G_clust.Edges.pairs(:,2);
    grainsA = grains.meanOrientation(grainIdA);
    grainsB = grains.meanOrientation(grainIdB);

    % extract options 
    sigma=opt.sigma;
    twin=opt.twin;
    nTwin=opt.nTwin;
    
    %Compute the max shear with principle stresses
    %Remember that principle stresses are eigenvalues and the ordering of
    %the principle stress is from large (1) to small (3). Look at mohrs 
    %circle for info.
    sigmaPrinciple = sigma.eig;
    tauMax = (max(sigmaPrinciple) - min(sigmaPrinciple)) / 2;  
    
    %Get array sizes 
    ngrains = length(G_clust.Nodes.Id);
    nedges = length(G_clust.Edges.pairs);
    
    %Allocate arrays for storage
    sigma13 = zeros(nedges,2);
    EffSFRelative = zeros(nedges,1);
    EffSF = zeros(nedges,2);
%     sym_ops = zeros(nedges,2);
    grainEffSF_tmp = zeros(nedges,2);
    typeEdge = G_clust.Edges.type;

    %Loop over edges
    for i = 1:nedges  
        %Only do calcs once for each edge
        if doCalc(i)
            type=typeEdge(i);
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
                grainEffSF_tmp(i,:) = [sigmaA(1,3) / (2*tauMax),sigmaB(1,3) / (2*tauMax)];
                EffSF(i,:) = [sigmaA(1,3) / (2*tauMax),sigmaB(1,3) / (2*tauMax)];

                EffSFRelative(i) = (sigmaA(1,3) - sigmaB(1,3)) / (2*tauMax);
            end
        end
    end %Loop over edges
    
    grainEffSF = zeros(ngrains, nTwin);
    for i = 1:nedges
        if doCalc(i)
            if typeEdge(i) <= nTwin && typeEdge(i) > 0
                grainEffSF(grainIdA(i),typeEdge(i))=grainEffSF_tmp(i,1);
                grainEffSF(grainIdB(i),typeEdge(i))=grainEffSF_tmp(i,2);
            end
        end
    end
    %Store arrays
    G_clust.Edges.sigma13(doCalc,:) = sigma13(doCalc,:);
    ind=unique(vertcat(grainIdA(doCalc),grainIdB(doCalc)));
    G_clust.Nodes.EffSF(ind,:) = grainEffSF(ind,:);
    G_clust.Edges.EffSFRelative(doCalc) = EffSFRelative(doCalc) ;
    G_clust.Edges.EffSF(doCalc,:) = EffSF(doCalc,:);
    G.Edges.sigma13(G.Edges.combineCleaned,:) = G_clust.Edges.sigma13;
    G.Nodes.EffSF(G_clust.Nodes.Id,:) = G_clust.Nodes.EffSF;
    G.Edges.EffSFRelative(G.Edges.combineCleaned) = G_clust.Edges.EffSFRelative;
    G.Edges.EffSF(G.Edges.combineCleaned,:) = G_clust.Edges.EffSF; 
    G.Edges.SFCalcd(G.Edges.combineCleaned) = true; 
    
    if opt.plot.do
        for i=1:nTwin
            %Plot 
            labelNodes=false;labelEdges=opt.plot.labelEdges; plotG=false; legendOn=opt.plot.legendOn;
            fhandle = plotGraph(grains,mergedGrains,G_clust,...
                G_clust.Nodes.EffSF(:,i),G_clust.Nodes.Id,...
                labelNodes,labelEdges,legendOn,plotG,[]);
            mtexColorbar;
            titleName=sprintf('EffSchmid, $%s$', twin{i}.name);
            mtexTitle(titleName);
        end
    end
    
end

