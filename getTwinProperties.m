function [twin] = getTwinProperties(twin)
%getTwinProperties Summary of this function goes here
    for i=1:length(twin)

        %Define the twin frame Rtw such that Rtw transforms crystal to twin
        twin{i}.Rtw=orientation.map(twin{i}.k1,twin{i}.CS.cAxis,...
            twin{i}.eta1,twin{i}.CS.aAxis); 
        
%         twin{i}.Rtw=getRtw(twin{i}.k1.hkl,twin{i}.eta1.UVTW,twin{i}.CS);
        %Define the transformation types
        tType{1}=orientation.byMatrix([-1  0  0;0 -1  0;0  0 1],twin{i}.CS);
        tType{2}=orientation.byMatrix([1  0  0;0 -1  0;0  0 -1],twin{i}.CS); 
        tType{3}=tType{1}*tType{2};
        tType{4}=tType{2}*tType{1};

        % Create the twin misorientation
        twin{i}.RMT=twin{i}.Rtw'*tType{twin{i}.actType}*twin{i}.Rtw;
        
        % Compute misorientation axis and angle
        twin{i}.axis=round(twin{i}.RMT.axis);
        twin{i}.angle=angle(twin{i}.RMT)/degree;
    end
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