function [Family] = GetFamily(ori,seg_angle)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    len=length(ori);
    oriSame=zeros(len,len,'logical');
    FamiliesComplete=zeros(len,1,'logical');
    Family=zeros(len,1);
    cnt=0;
    for i=1:length(ori)
        %Compute the misorientation between other grains
%         for j=1:length(ori)
        oriSame(i,:)=angle(ori(i),ori(:)) < seg_angle; %
%         end
        
        %check if oriSame(i,:) is part of any previous cluster 
        if i>1 
            for k=1:i-1
               if any((oriSame(i,:)+oriSame(k,:))==2)
                   oriSameCombined= logical(oriSame(i,:)+oriSame(k,:));
                   Family(oriSameCombined)=Family(k);
                   FamiliesComplete(oriSameCombined)=true;
               end
            end
        end
        
        if FamiliesComplete(oriSame(i,:))~=true
            cnt=cnt+1;
            FamiliesComplete(oriSame(i,:))=true;
            Family(oriSame(i,:))=cnt;
%             if cnt==9
%                fprintf('i=%d\n',i)
%             end
        end
        if sum(FamiliesComplete)==len
            break;
        elseif sum(FamiliesComplete)==len-1
            cnt=cnt+1;
            Family(~FamiliesComplete)=cnt;
            break;
        end 
    end
end

