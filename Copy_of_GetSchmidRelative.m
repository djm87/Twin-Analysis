function [sFRelative1,sFRelative2] = GetSchmidRelative(grains1,grains2,tt1,tt2,ct1,CRSS_TT1,CRSS_TT2,CRSS_CT1,sigma)
%The script returns the sign of the schmid factor for WE43 for given grains.
    CS=grains1.CS
    %% define slip systems and twinning systems
    %{10-12} extension twinning
    sSExtTwin1 = slipSystem(Miller(-1,0,1,1,CS,'uvtw'), Miller(1,0,-1,2,CS,'hkl'),CRSS_TT1);
    sSExtTwin1 = sSExtTwin1.symmetrise('antipodal')
  
    %{11-21} extension twinning
    sSExtTwin2 = slipSystem(Miller(1,1,-2,6,CS,'uvtw'), Miller(-1,-1,2,1,CS,'hkl'),CRSS_TT2);
    sSExtTwin2 = sSExtTwin2.symmetrise('antipodal');

    %{10-11} contraction twinning
    sSConTwin1 = slipSystem(Miller(2,-1,-1,-3,CS,'uvtw'), Miller(2,-1,-1,2,CS,'hkl'),CRSS_CT1);
    sSConTwin1 = sSConTwin1.symmetrise('antipodal');

    %Set the twin rotation axes for variant determination.
    twinAxes_tt1=tt1.axis.symmetrise; %Equivalent to using sSExtTwin1.n and 180 rotation
    twinAxes_tt2=tt2.axis.symmetrise;
    twinAxes_ct1=ct1.axis.symmetrise;
    %% Rotate loading to crystal
    rCS1=rotate(sigma,inv(grains1.meanOrientation))
    rCS2=rotate(sigma,inv(grains2.meanOrientation))
    
    for i=1:length(grains1)
       if grains1.prop.type(i)==1
        % Ext Twin 1 calculations
            for j=1:length(twinAxes_tt1)
                vari(j)=grains1(i).meanOrientation*twinAxes_tt1(j);

                %set a rotation around that axis
                rot(j)=rotation('axis',vari(j),'angle',tt1.angle);

                % rotate the c axis around the twin axis, and save orientation
                oriV(j)=rot(j)*grains1(i).meanOrientation;
                mis(j)=angle(oriV(j), grains2(i).meanOrientation)/degree;
            end
            [misIdeal(i),ActiveVar(i)]=min(abs(mis));
            
            sF1 = sSExtTwin1.SchmidFactor(rCS1(i));
            sF2 = sSExtTwin1.SchmidFactor(rCS2(i));
            % figure;plot(ebsd('Magnesium'),max(sFExtTwin1,[],2))

        elseif grains1.prop.type(i)==2
        % Ext Twin 2 calculations
            for j=1:length(twinAxes_tt2)
                vari(j)=grains1(i).meanOrientation*twinAxes_tt2(j);

                %set a rotation around that axis
                rot(j)=rotation('axis',vari(j),'angle',tt2.angle);

                % rotate the c axis around the twin axis, and save orientation
                oriV(j)=rot(j)*grains1(i).meanOrientation;
                mis(j)=angle(oriV(j), grains2(i).meanOrientation)/degree;
            end
            [misIdeal(i),ActiveVar(i)]=min(abs(mis));
            sF1 = sSExtTwin2.SchmidFactor(rCS1(i));
            sF2 = sSExtTwin2.SchmidFactor(rCS2(i));
            % figure;plot(ebsd('Magnesium'),max(sFExtTwin2,[],2))

        elseif grains1.prop.type(i)==3
        % Contraction Twin 1 calculations
            for j=1:length(twinAxes_tt2)
                vari(j)=grains1(i).meanOrientation*twinAxes_ct1(j);

                %set a rotation around that axis
                rot(j)=rotation('axis',vari(j),'angle',ct1.angle);

                % rotate the c axis around the twin axis, and save orientation
                oriV(j)=rot(j)*grains1(i).meanOrientation;
                mis(j)=angle(oriV(j), grains2(i).meanOrientation)/degree;
            end
            [misIdeal(i),ActiveVar(i)]=min(abs(mis));
            sF1 = sSConTwin1.SchmidFactor(rCS1(i));
            sF2 = sSConTwin1.SchmidFactor(rCS2(i));
            % figure;plot(ebsd('Magnesium'),max(sFConTwin1,[],2))
       elseif grains1.prop.type(i)==0
           sF1=1; sF2=0.5; ActiveVar(i)=1;
           
       end
        % compute the maximum relative Schmid factors
        sF1Active(i) = sF1(ActiveVar(i));      
        sF2Active(i) = sF2(ActiveVar(i));
    end
    
    for i=1:length(sF1Active)
        %Test to see if same max/min
        if  sF1Active(i)<0
            sFRelative1(i)=-1; %twin
        else
            sFRelative1(i)=1; %parent
        end
        
        if  sF2Active(i)<0
            sFRelative2(i)=-1; %twin
        else
            sFRelative2(i)=1; %parent
        end
        
        if  sFRelative2(i)==sFRelative1(i)
            sFRelative1(i)=0; %both negative or both positive
            sFRelative2(i)=0; 
        end
    end
    
end

