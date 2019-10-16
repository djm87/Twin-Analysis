function [group,typeout] = TestTwinRelationship(mori,tol,twin);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
    %Assume not twin and not to be grouped unless proven to be a twin type
    typeout=zeros(length(mori),1,'int8'); %non-zero describes twin type
    group=zeros(length(mori),1,'logical');
    cnt=0;
    for i=1:length(twin)
        for j=1:length(twin{i}.variantsToUse) %for double twins
            cnt=cnt+1;
            grouptmp=(angle(mori,twin{i}.RMT(j))' < tol)';
            typeout(grouptmp)=cnt; 
            group=logical(grouptmp+group); 
        end
    end
    
end


