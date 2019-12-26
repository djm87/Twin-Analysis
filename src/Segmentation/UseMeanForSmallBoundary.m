function [combine,GbInd] = UseMeanForSmallBoundary(pairs,GbInd,combine,grains,meanType,opt)
%UseMeanForSmallBoundary uses mean misorientation when less than GbMinLen
%boundary segments are present. The goal is subvert tripple junctions from
%causing unwanted grains merging. This feature has not been extensively
%tested.. 

%Edge segements are available.
    mineral=grains.mineral;
    gB=grains.boundary;
    gB_mineral = gB(mineral,mineral);
    gB_Id=gB_mineral.grainId;
    pairsCombine=pairs(combine,:);
    combineTmp=combine(combine);
    GbIndTmp=GbInd(combine);
    for i=1:length(pairsCombine)
        isTwinning=find(all(pairsCombine(i,:)==gB_Id,2)); %This is very expensive        
        if length(isTwinning)< opt.GbMinLen & meanType(i)==0
           combineTmp(i)=false;
        end
        GbIndTmp(i)=isTwinning(1);
    end
    GbInd(combine)=GbIndTmp;
    combine(combine)=combineTmp;

end