function [twinBoundary] = EBSDTwinBoundary(grains,mistol,nTwin,twin,mergeTripplePoints)
    
    mineral=grains.mineral;
    gB=grains.boundary;
    tripplePoints=unique(grains.triplePoints.boundaryId);
    if ~mergeTripplePoints
       gB(tripplePoints)=[];         
    end
    
    gB_mineral = gB(mineral,mineral);
    twinBoundary={};
    cnt=0;
    for i=1:nTwin 
        twinBoundary{i}=[];
        for j=1:length(twin{i}.variantsToUse) %for double twins
            cnt=cnt+1;
            isTwinning = angle(gB_mineral.misorientation,twin{i}.RMT(j)) < mistol(i);
            twinBoundary{cnt} = gB_mineral(isTwinning);
        end
    end
end