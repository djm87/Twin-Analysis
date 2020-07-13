function [Family] = GetFamily(ori,seg_angle)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%     len=length(ori);
%     oriSame=zeros(len,len,'logical');
%     FamiliesComplete=zeros(len,1,'logical');
%     Family=zeros(len,1);
%     cnt=1;
%     misOriAng={};
%     meanOri={};
    [Family,center] = calcCluster(ori,'maxAngle',seg_angle,'method','hierarchical');
    
%     figure; 
%     oR=fundamentalRegion(ori.CS)
%     plot(oR)
%     hold on 
%     plot(ori,ind2color(c))
%     caxis([1,5])
%     plot(center,'MarkerSize',10,'MarkerFaceColor','k','MarkerEdgeColor','k')
%     hold off
%     
%     while cnt<=len
%         %Compute the misorientation between other grains
%         ind=find(~FamiliesComplete,1,'first');
%         if isempty(ind)
%            break; 
%         end
%         %Get first orientation that isn't part of a family 
%         isOriSame=angle(ori(ind),ori(:)) < seg_angle; %  
%         
%         %Get mean of orientations that are within a tolerance
%         meanOri{cnt}=mean(ori(isOriSame));
%         
%         %Compute the misorientation
%         misOriAng{cnt}=angle(meanOri{cnt},ori(:));
%         isOriSame=misOriAng{cnt} < seg_angle; % 
%         if any(Family(isOriSame)~=0)
%             
%             uFamilies=unique(Family(isOriSame & Family~=0));
%             if length(uFamilies)>1
%                 warning('Too many families describe the families of similar tolerance.') ;
%                 warning('This could be a result of having a very large or too small Family seg angle');
%                 warning('Families will be merged');
%                 Family(isOriSame)=uFamilies(1); 
%                 FamiliesComplete(isOriSame)=true;
%                 for i=2:length(uFamilies)
%                     Family(uFamilies(i)==Family)=uFamilies(1);
%                 end
%             else
%                 Family(isOriSame)=uFamilies;
%                 FamiliesComplete(isOriSame)=true;
%             end
%         else
%             Family(isOriSame)=cnt;
%             FamiliesComplete(isOriSame)=true;
%             cnt=cnt+1;
%         end
%     end
%     
%     %Renumber the families
%     uFamilies=unique(Family);
%     FamilyTmp=Family;
%     for i=1:length(uFamilies)
%         Family(uFamilies(i)==FamilyTmp)=i;
%     end
    
%     for i=1:length(ori)
%         %Compute the misorientation between other grains
% %         for j=1:length(ori)
%         oriSame(i,:)=angle(ori(i),ori(:)) < seg_angle; %
% %         end
%         
%         %check if oriSame(i,:) is part of any previous cluster 
%         if i>1 
%             for k=1:i-1
%                if any((oriSame(i,:)+oriSame(k,:))==2)
%                    oriSameCombined= logical(oriSame(i,:)+oriSame(k,:));
%                    Family(oriSameCombined)=Family(k);
%                    FamiliesComplete(oriSameCombined)=true;
%                end
%             end
%         end
%         
%         if FamiliesComplete(oriSame(i,:))~=true
%             cnt=cnt+1;
%             FamiliesComplete(oriSame(i,:))=true;
%             Family(oriSame(i,:))=cnt;
% %             if cnt==9
% %                fprintf('i=%d\n',i)
% %             end
%         end
%         if sum(FamiliesComplete)==len
%             break;
%         elseif sum(FamiliesComplete)==len-1
%             cnt=cnt+1;
%             Family(~FamiliesComplete)=cnt;
%             break;
%         end 
%     end
end

