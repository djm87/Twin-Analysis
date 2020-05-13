function [mGrains] = CountTwins(G,grains,mGrains,opt)
%CountTwins produces twin count statistics 
    nGroups=max(G.Nodes.Group);
    twinCount = zeros(nGroups ,1);
    twinFamilyCount = zeros(nGroups ,1);
    for i=1:nGroups 
        ngroupId = find(i==G.Nodes.Group);
        nFamilyID=G.Nodes.FamilyID(ngroupId);
        nType = G.Nodes.type(ngroupId); 

        %Find the twins
        isTwin=nType>0& nType~=opt.twinUnknown;
        
        %If count the instances
        twinCount(i)=numel(nType(isTwin));
        twinFamilyCount(i)=numel(unique(nFamilyID(isTwin)));
    end

    mGrains.prop.twinCount=twinCount;
    mGrains.prop.twinFamilyCount=twinFamilyCount;
end

