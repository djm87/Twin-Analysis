function [group,type] = TestTwinRelationship(mori,tol,twin,type);

    %Assume not twin and not to be grouped unless proven to be a twin type
    group=zeros(length(mori),1,'logical');
    cnt=0;
    for i=1:length(twin)-1
        for j=1:length(twin{i}.variantsToUse) %for double twins
            cnt=cnt+1;
            grouptmp=(angle(mori,twin{i}.RMT(j))' < tol)';
            type(grouptmp)=cnt; 
            group=logical(grouptmp+group); 
        end
    end
    group(type==length(twin))=true; %auto group unknown twin type
end


