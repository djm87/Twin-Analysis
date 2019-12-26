function [group,type] = TestTwinRelationship(mori,tol,opt,type);

    %Assume not twin and not to be grouped unless proven to be a twin type
    group=zeros(length(mori),1,'logical');
    cnt=0;
    for i=1:opt.nTwin
        for j=1:length(opt.twin{i}.variantsToUse) %for double twins
            cnt=cnt+1;
            grouptmp=(angle(mori,opt.twin{i}.RMT(j))' < tol(i))';
            type(grouptmp)=cnt; 
            group=logical(grouptmp+group); 
        end
    end
    group(type==opt.twinUnknown)=true; %auto group unknown twin type
end


