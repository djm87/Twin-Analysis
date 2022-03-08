function [Family] = GetFamily(ori,misTol,LABTol)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%     len=length(ori);
%     oriSame=zeros(len,len,'logical');
%     FamiliesComplete=zeros(len,1,'logical');
%     Family=zeros(len,1);
%     cnt=1;
%     misOriAng={};
%     meanOri={};
    [cId,fCenters] = calcCluster(ori,'maxAngle',misTol,'method','hierarchical');

    %Get the misorientation between family
    while 1
        FamilyMis=round(triu(angle_outer(fCenters,fCenters)./degree),4);
        [r,c]=find(FamilyMis~=0);
 
        if isempty(r) || isempty(c), break; end
        unchangedCnt=0;
        for j=1:size(r,1)
            rid=find(cId==r(j));
            cid=find(cId==c(j));
            if ~isempty(rid) && ~isempty(cid)
                orir=ori(rid);
                oric=ori(cid);

                d = full(abs(dot_outer(orir,oric)));
                dr = full(abs(dot_outer(orir,orir)));
                dc = full(abs(dot_outer(oric,oric)));

                dstart=d;
                cIdstart=cId;
                % progress(0,length(ori));
                while 1

                  % find smallest pair
                  [omega,id] = max(d(:));

                  if omega<cos(LABTol), break; end

                  [l,m] = ind2sub(size(d),id);

                  cId(cid(m))=cId(rid(l));

                  d(l,m)=0;
                end
                if all(d(:)==dstart(:)) && all(cId(:)==cIdstart(:))
                    unchangedCnt=unchangedCnt+1;
                end

            else
                unchangedCnt=unchangedCnt+1;
            end
        end
        if unchangedCnt==size(r,1), break; end

%         fcnt=1;
%         for j=1:max(cId)
%             ind=cId==j;
%             if any(ind)
%                 fCenters(fcnt)=mean(ori(ind));
%                 fcnt=fcnt+1;
%             end
%         end
    end
    Family=cId;
    
    cnt=1;
    for i=1:max(Family)
        ind=i==Family;
        if any(ind)
            Family(ind)=cnt;
            cnt=cnt+1;
        end
    end
%     fCenters(fcnt:end)=[];
    
    
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

