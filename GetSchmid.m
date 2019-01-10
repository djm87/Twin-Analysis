function [sFRelativeMaxMag] = GetSchmid(grains,CRSS_Basal,CRSS_Pris,CRSS_pyr,CRSS_TT1,CRSS_TT2,CRSS_CT1,loadDir)
%The script returns the sign of the schmid factor for WE43 for given grains.
    CS=grains.CS
    %% define slip systems and twinnings
    sSBasal = slipSystem(Miller(2,-1,-1,0,CS,'uvtw'), Miller(0,0,0,1,CS,'hkl'),CRSS_Basal/CRSS_Basal);
    sSBasal = sSBasal.symmetrise%('antipodal')

    sSPrismatic = slipSystem(Miller(-1,2,-1,0,CS,'uvtw'), Miller(1,0,-1,1,CS,'hkl'),CRSS_Pris/CRSS_Basal);
    sSPrismatic = sSPrismatic.symmetrise%('antipodal')

    % type II pyramidal
    sSPyramidal = slipSystem(Miller(2,-1,-1,-3,CS,'uvtw'), Miller(2,-1,-1,2,CS,'hkl'),CRSS_pyr/CRSS_Basal);
    sSPyramidal = sSPyramidal.symmetrise%('antipodal')

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
    if strcmp(class(grains),'EBSD')
        rCS = grains.orientations \ loadDir;
    else
        rCS = grains.meanOrientation \ loadDir;
    end
    
    %% Basal slip calculations

    sFBasal = abs(sSBasal.SchmidFactor(rCS));

    % figure;plot(ebsd('Magnesium'),max(sFBasal,[],2))

    %% Prismatic slip calculations

    sFPrismatic = abs(sSPrismatic.SchmidFactor(rCS));
    % figure;plot(ebsd('Magnesium'),max(sFPrismatic,[],2))

    %% Pyramidal slip calculations

    sFPyramidal = abs(sSPyramidal.SchmidFactor(rCS));
    % figure;plot(ebsd('Magnesium'),max(sFPyramidal,[],2))

    %% Ext Twin 1 calculations

    sFExtTwin1 = abs(sSExtTwin1.SchmidFactor(rCS));
    % figure;plot(ebsd('Magnesium'),max(sFExtTwin1,[],2))

    %% Ext Twin 2 calculations

    sFExtTwin2 = abs(sSExtTwin2.SchmidFactor(rCS));
    % figure;plot(ebsd('Magnesium'),max(sFExtTwin2,[],2))

    %% Contraction Twin 1 calculations

    sFConTwin1 = abs(sSConTwin1.SchmidFactor(rCS));
    % figure;plot(ebsd('Magnesium'),max(sFConTwin1,[],2))

    %% compare SchmidFactors scaled by CRSS

    % combine all slip systems
%     sS = [sSBasal;sSPrismatic;sSPyramidal;sSExtTwin1;sSExtTwin2;sSConTwin1];
    sS = [sSExtTwin1];
    sFRelative=sS.SchmidFactor(rCS);
    % compute the relative Schmid factors scaled by CRSS
%     sFRelative = sS.SchmidFactor(rCS,'relative');

    %% compute the maximum relative Schmid factors
    [sFmax,sFidmax] = max(sFRelative,[],2);    
    [sFmin,sFidmin] = min(sFRelative,[],2);

    for i=1:length(sFmax)
        %Test to see if same max/min
        if sFmax(i)+sFmin(i) >0
            sFRelativeMaxMag(i)=1; 
        elseif sFmax(i)+sFmin(i) <0
            sFRelativeMaxMag(i)=-1; 
        else 
            sFRelativeMaxMag(i)=1; %Could go either way
        end

    end
    
%     [sFmax,sFidmax]=max(abs(sFRelative),[],2); 
%     for i=1:length(sFmax)
%         sFRelativeMaxMag(i)=sFRelative(i,sFidmax(i));
%     end
end

