function [G_Family,G_clust] = GetSchmidRelative(G_Family,G_clust,oriFamily,typeEdge,grains,mGrains,opt)
%The script returns the effective schmid factor on K1 plane in the eta1
%direction. In addition, the symmetry operators are returned. Since this is
%called before we know the parent, we don't know the sign of the twin
%variant rotation axis. For this reason variant extraction is done in a
%seperate routine after the family tree is formed.
    
    %Get array sizes 
    nedges = size(G_Family.Edges.familyPair,1);
    
    % extract mean orientation for each family edge relationship in group
    grainsA = oriFamily(:,1);
    grainsB = oriFamily(:,2);

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
    
    %Allocate arrays for storage
    sigma13 = zeros(nedges,2);
    EffSF = zeros(nedges,2);

    %Loop over edges
    for i = 1:nedges  
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
            EffSF(i,:) = [sigmaA(1,3) / (2*tauMax),sigmaB(1,3) / (2*tauMax)];
        end
    end %end loop over edges
    
    %Store arrays
    G_Family.Edges.sigma13 = sigma13;
    G_Family.Edges.EffSF = EffSF;
    
    %Store results in global graph for plotting
    G_clust.Nodes.EffSF=zeros(length(G_clust.Nodes.Group),nTwin);
    for i=1:length(mGrains) 
        egroupFId = find(i==G_Family.Edges.Group);
        if ~isempty(egroupFId)
            ngroupId = find(i==G_clust.Nodes.Group);
            nFamily = G_clust.Nodes.FamilyID(ngroupId);     
            familyPair=G_Family.Edges.familyPair(egroupFId,:);
            EffSFtmp=EffSF(egroupFId,:);
            type=typeEdge(egroupFId);

            nEffSF=zeros(length(ngroupId),nTwin);
            nEffSFType=zeros(length(ngroupId),1);
            
            for k=1:opt.nTwin
                %find relationships of type k
                [row]=find(type==k);
                familyPairsub=familyPair(row,:);
                EffSFsub=EffSFtmp(row,:);
                [famList]=unique(familyPairsub);
                for j=1:length(famList)
                    nEffSF(nFamily==famList(j),k)=mean(EffSFsub(familyPairsub==famList(j)));
                end
                G_clust.Nodes.EffSF(ngroupId,k)=nEffSF(:,k);
            end            
        end
    end
        
    if opt.plot.do
        for i=1:nTwin
            %Plot 
            labelNodes=false;labelEdges=opt.plot.labelEdges; plotG=false; legendOn=opt.plot.legendOn;
            fhandle = plotGraph(grains,mGrains,G_clust,...
                G_clust.Nodes.EffSF(:,i),G_clust.Nodes.Id,...
                labelNodes,labelEdges,legendOn,plotG,[]);
            mtexColorbar;
            titleName=sprintf('EffSchmid, $%s$', twin{i}.name);
            mtexTitle(titleName);
        end
    end
    
end

