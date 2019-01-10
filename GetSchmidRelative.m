function G = GetSchmidRelative(G,twins,CRSS,sS,sigma)
%The script returns the sign of the schmid factor for WE43 for given grains.
    % extract grains 
    grainId1=G.Edges.pairs(:,1);
    grainId2=G.Edges.pairs(:,2);
    grains1=G.Nodes.meanOrientation(grainId1);
    grains2=G.Nodes.meanOrientation(grainId2);
    
    % Rotate loading to crystal frame  
    rCS1=rotate(sigma,inv(grains1));
    rCS2=rotate(sigma,inv(grains2));
             
    %Set the twin rotation axes for variant determination.
    twinAxes={};
    for i=1:length(twins)
        twinAxes{i}=twins{i}.axis.symmetrise;
    end
    
    %Find twin shear direction 
    
    
    %Initialize variables
    nedges=length(G.Edges.pairs); 
    TwinVarActive=zeros(nedges,1);
    sF1Active=zeros(nedges,1);
    sF2Active=zeros(nedges,1);    
    sFRelative12=zeros(nedges,1); 
    sFActiveVar=zeros(nedges,1); 
    for i=1:nedges 
        %For some pair (edge) get the type of twin relation
        type=G.Edges.type(i); %vectorize this
        
        if type >0 %0 if not twin type
            %For the first grain in the pair compute the twin axes
            vari=grains1(i)*twinAxes{type};

            %set a rotation around that axis
            rot=rotation('axis',vari,'angle',twins{type}.angle);

            % rotate around the twin axis, and save orientation
            oriV=rot*grains1(i);
            
            %calculate misorientation between twin variants in the second
            %grain
            mis=angle(oriV, grains2(i))/degree; %all the varients of the twin mode

            [misIdeal(i),TwinVarActive(i)]=min(abs(mis)); %active twin variant
            
            sF1 = sS{type}.SchmidFactor(rCS1(i)); %Compute Schmid for the two twin orientations
            sF2 = sS{type}.SchmidFactor(rCS2(i));
            
            %Look at: https://github.com/mtex-toolbox/mtex/issues/288 and 
            %https://gist.github.com/jhiscocks/4acc046799d90dfe54e01936527a51e5
            %to determine the appropriate general way of calculating this. 
            %also see rod's article. a given twin should have one slip
            %direction and plane.
%             sS{type}.SchmidFactor(rotate(sigma,inv(oriV(TwinVarActive(i)))))
            
            % compute the maximum relative Schmid factors
            sF1Active(i) = sF1(TwinVarActive(i));      
            sF2Active(i) = sF2(TwinVarActive(i));
            sFRelative12(i)=sF1Active(i)-sF2Active(i);
        end
    end
    
    %Save to the graph data structure
    G.Edges.SF=[sF1Active,sF2Active]; %Schmid factors for relation type 
    G.Edges.SFRelative12=sFRelative12;
    G.Edges.SFActiveVar=sFActiveVar; 
    
end

