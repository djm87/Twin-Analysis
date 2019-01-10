function G = GetSchmidRelative(G,twins,CRSS,sS,sigma)
%The script returns the sign of the schmid factor for WE43 for given grains.
    % extract grains 
    grainId1=G.Edges.pairs(:,1);
    grainId2=G.Edges.pairs(:,2);
    grains1=G.Nodes.meanOrientation(grainId1);
    grains2=G.Nodes.meanOrientation(grainId2);
    
    %Set the twin rotation axes for variant determination.
    twinAxes={};
    for i=1:length(twins)
        twinAxes{i}=twins{i}.axis.symmetrise;
    end
    
    % Rotate loading to crystal frame  
    rCS1=rotate(sigma,inv(grains1));
    rCS2=rotate(sigma,inv(grains2));
    
    %Initialize variables
    ngrains=length(grains1); 
    sF1Active=zeros(ngrains,1);
    sF2Active=zeros(ngrains,1);    
    sFRelative12=zeros(ngrains,1); 
    sFActiveVar=zeros(ngrains,1); 
    for i=1:ngrains
        
        type=G.Edges.type(i); %vectorize this
       if type>0 %i.e. twin relation
%            tic
            vari=grains1(i)*twinAxes{type};

            %set a rotation around that axis
            rot=rotation('axis',vari,'angle',twins{type}.angle);

            % rotate the c axis around the twin axis, and save orientation
            oriV=rot*grains1(i);
            mis=angle(oriV, grains2(i))/degree;
%            toc
            [misIdeal(i),sFActiveVar(i)]=min(abs(mis));
            
            sF1 = sS{type}.SchmidFactor(rCS1(i));
            sF2 = sS{type}.SchmidFactor(rCS2(i));

            % compute the maximum relative Schmid factors
            sF1Active(i) = sF1(sFActiveVar(i));      
            sF2Active(i) = sF2(sFActiveVar(i));
            sFRelative12(i)=sF1Active(i)-sF2Active(i);
       end
    end
    
    %Save to the graph data structure
    G.Edges.SF=[sF1Active,sF2Active]; %Schmid factors for relation type 
    G.Edges.SFRelative12=sFRelative12;
    G.Edges.SFActiveVar=sFActiveVar; 
    
end

