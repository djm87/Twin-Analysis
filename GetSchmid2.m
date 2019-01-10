function [sFRelativeMaxMag] = GetSchmid(grains,CRSS_Basal,CRSS_Pris,CRSS_pyr,CRSS_TT1,CRSS_TT2,CRSS_CT1,sigma,type)
%The script returns the sign of the schmid factor for WE43 for given grains.
    CS=grains.CS
    %% define slip systems and twinning systems
    %{10-12} extension twinning
    sSExtTwin1 = slipSystem(Miller(-1,0,1,1,CS,'uvtw'), Miller(1,0,-1,2,CS,'hkl'),CRSS_TT1/CRSS_Basal);
    sSExtTwin1 = sSExtTwin1.symmetrise%('antipodal')

    %{11-21} extension twinning
    sSExtTwin2 = slipSystem(Miller(1,1,-2,6,CS,'uvtw'), Miller(-1,-1,2,1,CS,'hkl'),CRSS_TT2/CRSS_Basal);
    sSExtTwin2 = sSExtTwin2.symmetrise%('antipodal');

    %{10-11} contraction twinning
    sSConTwin1 = slipSystem(Miller(2,-1,-1,-3,CS,'uvtw'), Miller(2,-1,-1,2,CS,'hkl'),CRSS_CT1/CRSS_Basal);
    sSConTwin1 = sSConTwin1.symmetrise%('antipodal');

    %% Rotate loading to crystal
    rCS=rotate(sigma,inv(grains.meanOrientation))
    
    if grains.prop.type==1
    % Ext Twin 1 calculations
        sF = sSExtTwin1.SchmidFactor(rCS);
        % figure;plot(ebsd('Magnesium'),max(sFExtTwin1,[],2))

    % Ext Twin 2 calculations
    elseif grains.prop.type==2
        sF = sSExtTwin2.SchmidFactor(rCS);
        % figure;plot(ebsd('Magnesium'),max(sFExtTwin2,[],2))

    % Contraction Twin 1 calculations
    elseif grains.prop.type==3
        sF = sSConTwin1.SchmidFactor(rCS);
        % figure;plot(ebsd('Magnesium'),max(sFConTwin1,[],2))
    end
    %% compare SchmidFactors scaled by CRSS

    %% compute the maximum relative Schmid factors
    [sFmax,sFidmax] = max(sF,[],2);    
    [sFmin,sFidmin] = min(sF,[],2);
    varient=
    for i=1:length(sFmax)
        %Test to see if same max/min
        if sFmax(i)+sFmin(i) >0
            sFRelative(i)=1; 
        elseif sFmax(i)+sFmin(i) <0
            sFRelative(i)=-1; 
        else 
            sFRelative(i)=0; %Not informative
        end
    end
    
end

