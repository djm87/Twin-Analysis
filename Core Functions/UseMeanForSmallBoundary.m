function combine = UseMeanForSmallBoundary(pairs,combine,grains,typeMeanMisTol,minNEdgeMistol)
%UseMeanForSmallBoundary uses mean misorientation when less than minNEdgeMistol
%Edge segements are available.
    mineral=grains.mineral;
    gB=grains.boundary;
    gB_mineral = gB(mineral,mineral);
    gB_Id=gB_mineral.grainId;
    pairsCombine=pairs(combine,:);
    combineTmp=combine(combine);
    parfor i=1:length(pairsCombine)   
        isTwinning=find(all(pairsCombine(i,:)==gB_Id,2));
        if length(isTwinning)<minNEdgeMistol & typeMeanMisTol(i)==0
           combineTmp(i)=false;
        end
    end
    combine(combine)=combineTmp;

end