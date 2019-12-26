function [twin] = getTwinProperties(twin)
%getTwinProperties Summary of this function goes here
    
    ntwin=length(twin);
    %Add twin for unidentified internal grains 
    twin{ntwin+1}.name='unknown internal';
    ntwin=ntwin+1;
    toRemove=zeros(ntwin,1,'logical');
    newTwinCnt=1;
    for i=1:ntwin-1
        
        %Define the transformation types
        tType{1}=orientation.byMatrix([-1  0  0;0 -1  0;0  0 1],twin{i}.CS);
        tType{2}=orientation.byMatrix([1  0  0;0 -1  0;0  0 -1],twin{i}.CS);   
        
        %handle double and single twins
        if length(twin{i}.k1)==2
           %Define the twin frame Rtw such that Rtw transforms crystal to twin
            twin{i}.Rtw(1)=orientation.map(twin{i}.k1(1),twin{i}.CS.cAxis,...
                twin{i}.eta1(1),twin{i}.CS.aAxis); 
            twin{i}.Rtw(2)=orientation.map(twin{i}.k1(2),twin{i}.CS.cAxis,...
                twin{i}.eta1(2),twin{i}.CS.aAxis);
            
            % Create the twin misorientation
            twin{i}.RMT=twin{i}.Rtw(1)'*tType{twin{i}.actType}*twin{i}.Rtw(1)*twin{i}.CS*...
                twin{i}.Rtw(2)'*tType{twin{i}.actType}*twin{i}.Rtw(2);
            
            %Extract
            [uAngle,IA,IC]=unique(round(twin{i}.RMT.angle,2));
            twin{i}.RMT=twin{i}.RMT(IA(twin{i}.variantsToUse));
            
            % Compute misorientation axis and angle
            twin{i}.axis=round(twin{i}.RMT.axis);
            twin{i}.angle=angle(twin{i}.RMT)/degree;
            
            twin{ntwin+newTwinCnt}=twin{i};
            twin{ntwin+newTwinCnt}.RMT=twin{i}.RMT;
            twin{ntwin+newTwinCnt}.angle=twin{i}.angle;
            twin{ntwin+newTwinCnt}.axis=twin{i}.axis;
            twin{ntwin+newTwinCnt}.sS=slipSystem(twin{i}.eta1(1),twin{i}.k1(1));
            twin{ntwin+newTwinCnt}.axisVariants=twin{ntwin+1}.axis.symmetrise;
            twin{ntwin+newTwinCnt}.name=twin{i}.name;
            twin{ntwin+newTwinCnt}.variantsToUse=1;
            newTwinCnt=newTwinCnt+1;
            toRemove(i)=true;
        else
            %Define the twin frame Rtw such that Rtw transforms crystal to twin
            twin{i}.Rtw=orientation.map(twin{i}.k1,twin{i}.CS.cAxis,...
                twin{i}.eta1,twin{i}.CS.aAxis); 
            
            % Create the twin misorientation
            twin{i}.RMT=twin{i}.Rtw'*tType{twin{i}.actType}*twin{i}.Rtw;

            % Compute misorientation axis and angle
            twin{i}.axis=round(twin{i}.RMT.axis);
            twin{i}.angle=angle(twin{i}.RMT)/degree;
            
            % Compute twin axis variants 
            twin{i}.axisVariants=twin{i}.axis.symmetrise;
            
            % Compute the slip systems for schmid calculations
            twin{i}.sS=slipSystem(twin{i}.eta1,twin{i}.k1);

        end
    end
    twin(find(toRemove))=[];
end

% function [Rtw] = getRtw(k1in,eta1in,CS)

%     a = 3.23;  alpha = CS.alpha; 
%     b = 3.23;  gamma = CS.gamma;  
%     c = 5.15;  beta = CS.beta;
%     
%     av = [a 0 0]; bv = [b*cos(gamma) b*sin(gamma) 0];
%     cv= [c*cos(beta) c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma) ...
%         c*sqrt(1+2*cos(alpha)*cos(beta)*cos(gamma)-cos(alpha)^2-...
%         cos(beta)^2-cos(gamma)^2)/sin(gamma)];
%     L=CS.axes.xyz';
%     av=L(:,1)';bv=L(:,2)';cv=L(:,3)';
%     
% %     L=[av' bv' cv'];
%     Q=[cross(bv,cv)' cross(cv,av)' cross(av, bv)']./det(L);
%     
%     %eta1 twin frame
%     
%     eta1=[eta1in(1)-eta1in(3); eta1in(2)-eta1in(3); eta1in(4)];
%     eta1_twin=L*eta1;
%     eta1_twin=eta1_twin./norm(eta1_twin,2);
%     
%     %k1 twin frame
%     k1=k1in([1,2,4])';
%     k1_twin=Q*k1;
%     k1_twin=k1_twin./norm(k1_twin,2);
%     s_twin=cross(k1_twin,eta1_twin);
%     
%     %Fill the rows of Rtw and convert to orientation
%     RtwMatrix=[eta1_twin';s_twin';k1_twin'];
%     Rtw=orientation.byMatrix(RtwMatrix);
% 
% end