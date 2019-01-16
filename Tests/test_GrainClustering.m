for i=246%length(parentId)
    grainCluster=grains(parentId==mergedGrains(i).id); 
    
    %group similar orientations
    len=length(grainCluster);
    oriSame=zeros(len,len,'logical');
    FamiliesComplete=zeros(len,1,'logical');
    cnt=1;
    for j=1:len
        for k=1:len
            oriSame(j,k)=angle(grainCluster(j).meanOrientation,grainCluster(k).meanOrientation) /degree <10
        end
        
        if FamiliesComplete(oriSame(j,:))~=true
            FamiliesComplete(oriSame(j,:))=true
            grains.prop.familyID(grainCluster.id(oriSame(j,:)))=cnt;
            cnt=cnt+1;
        end
        if sum(FamiliesComplete)==len
            break;
        elseif sum(FamiliesComplete)==len-1
            grains.prop.familyID(grainCluster.id(~FamiliesComplete))=cnt;
            break;
        end
    end
     grains.prop.familyID(grainCluster.id)
%     grainCluster.id(iSame
%     grains.prop.familyID=zeros(length(grains),1,'logical');
%     
%     grains.prop.type(i)=tmpType; 
    %grains(grainCluster.id) returns the same thing as
    %grains(parentId==mergedGrains(i).id) so we can write
    %family result to global grains.
    
    
%     figure
%     plot(ebsd_merged(mergedGrains),ebsd_merged(mergedGrains).orientations)    
%     hold on
%     plot(mergedGrains.boundary,'linewidth',2)
%     hold off 
    
%     for j=1:length(mergedGrains(i))
%         if length(mergedGrains(i))>1
%            'test'
%         end
%     end
end