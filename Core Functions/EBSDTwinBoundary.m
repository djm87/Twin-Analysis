function [twinBoundary] = EBSDTwinBoundary(grains,Mistol,twin)
    ntwins=length(twin);
    mineral=grains.mineral;
    gB=grains.boundary;
    gB_mineral = gB(mineral,mineral);
    twinBoundary={};
    cnt=0;
    for i=1:ntwins-1 
        twinBoundary{i}=[];
        for j=1:length(twin{i}.variantsToUse) %for double twins
            cnt=cnt+1;
            isTwinning = angle(gB_mineral.misorientation,twin{i}.RMT(j)) < Mistol;
            twinBoundary{cnt} = gB_mineral(isTwinning);
        end
    end
end