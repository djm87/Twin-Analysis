function [group,type] = TestTwinRelationship(mori,tol,opt,type)

    %Assume not twin and not to be grouped unless proven to be a twin type
    group=zeros(length(mori),1,'logical');
    cnt=0;
    for i=1:opt.nMori
        for j=1:length(opt.mori{i}.variantsToUse) %for double twins
            cnt=cnt+1;
            grouptmp=(angle(mori,opt.mori{i}.RMT(j))' < tol(i))';
            type(grouptmp)=cnt; 
            group=logical(grouptmp+group); 
        end
    end
    group(type==opt.moriUnknown)=true; %auto group unknown twin type
end


