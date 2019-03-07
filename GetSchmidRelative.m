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
    isAxisVariant = zeros(nedges,1,'uint8');
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
%                 tic
%                 angle(ori1, ori2.symmetrise,'noSymmetry')./ degree
                
                mori=inv(grainsA(i).symmetrise) *grainsB(i).symmetrise;
                MM=mori(:,:)*twin{type}.RMT'; 
%                 for j=1:144
%                     angs(j)=acos((trace(MM(j).matrix)-1)/2)/degree
%                 end
                angs=angle(MM, 'noSymmetry') ./ degree;
                [~,id]=min(angs);
                [sym_ops(i,1),sym_ops(i,2)] = ...
                    ind2sub([length(CS) length(CS)],id);
                
                %Determine the active twin variant
                actTwinAxis=CS(sym_ops(i,1))*twin{type}.RMT.axis;
                isVar=vector3d(round(actTwinAxis))==vector3d(twin{type}.axisVariants);
                assert(and(~isempty(isVar),sum(isVar)==1),'Twin variant not in list!');
                isAxisVariant(i)=find(isVar);
%                 L=CS.matrix
%                 Lg2=inv(L.*grainsA(i).matrix);
%                 g1g2tLt=(Lg2*grainsB(i));
%                 M = CS*g1g2tLt;
% %                 Find the difference between twin and grain misorientation
%                 MM = twin{type}.RMT'*M;
% %                 toc
%                 angs2 = angle(MM, 'noSymmetry')' / degree;
%                 [angs',angs2]
%                 %Extract the closest match along with the symmetry operations
%                 [~,id] = min(angs2);
%                 [sym_ops(i,1),sym_ops(i,2)] = ...
%                     ind2sub([length(CS) length(CS)],id);
%                 
%                 actTwinAxis=CS(sym_ops(i,1))*twin{type}.RMT.axis

                %Find the resolved shear stress on the K1 plane in the eta1
                %direction
%                 if length(twin{type}.Rtw)==2
%                     aA=twin{type}.Rtw(1).matrix * CS(sym_ops(i,1)).matrix * grainsA(i).matrix;
%                     aB=twin{type}.Rtw(2).matrix * CS(sym_ops(i,2)).matrix * grainsB(i).matrix;
%                 else
%                     aA=twin{type}.Rtw.matrix * CS(sym_ops(i,1)).matrix * inv(grainsA(i).matrix);
%                     aB=twin{type}.Rtw.matrix * CS(sym_ops(i,2)).matrix * inv(grainsB(i).matrix);  
%                 end
%                 sigma13(i,1) = dot(aA(1,:), aA(3,:) * sigma.matrix');
%                 sigma13(i,2) = dot(aB(1,:), aB(3,:) * sigma.matrix');
                
                gA=inv(grainsA(i).symmetrise); 
                gB=inv(grainsB(i).symmetrise);

                sigmaA = matrix(twin{type}.Rtw * (gA(sym_ops(i,1)) * sigma));
                sigmaB = matrix(twin{type}.Rtw * (gB(sym_ops(i,2)) * sigma));
                sigma13(i,1)=sigmaA(1,3);
                sigma13(i,2)=sigmaB(1,3);
                
                grainEffSF(grainIdA(i),type) = sigma13(i,1) / (2*tauMax);
                grainEffSF(grainIdB(i),type) = sigma13(i,2) / (2*tauMax);
                               
%                 sSTwin = slipSystem(twin{type}.k1,twin{type}.eta1).symmetrise
%                 RSCA= grainsA(i) \ sigma;
%                 RSCB=grainsB(i) \ sigma;
%                 schmid=sSTwin.SchmidFactor(RSCA)   
%                 sSTwin.SchmidFactor(RSCB)                
%                 
%                 Lg2=CS(sym_ops(i,2))*grainsB(i);
%                 g1g2tLt=(grainsA(i)*Lg2');
% 
%                 M = CS(sym_ops(i,1)) * g1g2tLt;
%                 round(M.axis('nosymmetry'))
%                 M.angle/degree
%                 var=CS(sym_ops(i,1))*twin{type}.RMT*CS(sym_ops(i,2))
%                 
% %                 test.angle('nosymmetry')/degree
% %                 test2=round(test.axis('nosymmetry'))
% %                 [test2.h test2.k test2.i  test2.l reshape(test.angle('nosymmetry')/degree,[144,1])] 
% 
%                 
                %Verify twin variant selection 
%                 twinAxes=twin{type}.axis.symmetrise;
%                 vari=grainsA(i)*twinAxes;
%                 %set a rotation around that axis
%                 rot=rotation('axis',vari,'angle',twin{type}.angle);
% 
%                 % rotate around the twin axis,
%                 oriV=rot*grainsA(i);
%                 %calculate misorientation between twin variants in the second
%                 %grain
%                 mis=angle(oriV, grainsB(i))/degree; %all the varients of the twin mode
% 
%                 [misIdeal,TwinVarActive]=min(abs(mis)); %active twin variant
%                 twinAxes(TwinVarActive)
%                 tic
            end
        end %Loop over edges
    end %Loop over twin types
    
    %Store arrays
    G.Nodes.EffSF = grainEffSF;
    G.Edges.sigma13 = sigma13;
    G.Edges.isAxisVariant=isAxisVariant;
    G.Edges.EffSF = sigma13 ./ (2 * tauMax);
    G.Edges.EffSFRelative = G.Edges.EffSF(:,1) - G.Edges.EffSF(:,2);
    
end

