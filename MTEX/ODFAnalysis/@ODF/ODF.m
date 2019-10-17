classdef ODF < dynOption
  
  properties    
    components = {};     % the ODF components, e.g., unimodal, fibres,
    weights    = [];     % the weights
  end
  
  properties (Dependent = true)
    CS          % crystal symmetry for ODF
    SS          % specimen symmetry for ODF
    antipodal   % mori =? inv(mori)
    bandwidth   % harmonic degree
  end
  
  methods
    function odf = ODF(components,weights)
      
      if nargin == 0, return;end
          
      odf.components = ensurecell(components);
      if nargin == 2
        odf.weights    = weights;
      else
        odf.weights = 1;
      end      
    end
    
    function odf = set.CS(odf,CS)
      for i = 1:numel(odf.components)
        odf.components{i}.CS = CS;
      end
    end
    
    function CS = get.CS(odf)
      if ~isempty(odf.components)
        CS = odf.components{1}.CS;
      else
        CS = [];
      end      
    end
    
    function odf = set.SS(odf,SS)
      for i = 1:numel(odf.components)
        odf.components{i}.SS = SS;
      end
    end
    
    function SS = get.SS(odf)
      if ~isempty(odf.components)
        SS = odf.components{1}.SS;
      else
        SS = [];
      end      
    end
    
    function odf = set.antipodal(odf,antipodal)
      for i = 1:numel(odf.components)
        odf.components{i}.antipodal = antipodal;
      end
    end
    
    function antipodal = get.antipodal(odf)
      if ~isempty(odf.components)
        antipodal = odf.components{1}.antipodal;
      else
        antipodal = false;
      end      
    end
    
    function L = get.bandwidth(odf)
      L = 0;
      for i = 1:numel(odf.components)
        L = max(L,odf.components{i}.bandwidth);
      end
    end
    
    function odf = set.bandwidth(odf,L)
      for i = 1:numel(odf.components)
        odf.components{i}.bandwidth =L;
      end
    end
    
  end
  
  methods (Static = true)
    
    function [odf,resvec] = interp(ori,values,varargin)
      % compute an ODF by interpolating orientations and weights

      % construct the uniform portion first
      values = values(:);
      m = min(values);
      values = values - m;
      
      odf = m * uniformODF(ori.CS,ori.SS);
      
      % grid for representing the ODF
      res = get_option(varargin,'resolution',3*degree);
      exact = get_option(varargin,'exact',false);
      if exact
          S3G=ori
      else
          S3G = equispacedSO3Grid(ori.CS,ori.SS,'resolution',res);
      end
%       S3G=ori
%       S3G = regularSO3Grid(ori.CS,ori.SS,'resolution',10*degree)     % specify the resolution
      % kernel for representing the ODF
      psi = get_option(varargin,'kernel',deLaValeePoussinKernel('halfwidth',res));

      % system matrix
      M = psi.K_symmetrised(S3G,ori,ori.CS,ori.SS);
      
%       x = lsqlin(C,d,A,b,[],[],[],[],[],options)
%       options = optimoptions('lsqlin','Algorithm','interior-point',...
%           'Display','iter','ConstraintTolerance',1e-8,'OptimalityTolerance',1e-8);
%       [w,resnorm,residual,exitflag,output,lambda] = lsqlin(M',values,[],[],[],[],[],[],[],options);
      
%       options = optimoptions('TolX',1e-3);
%       [w,resnorm,residual,exitflag,output,lambda] = lsqnonneg(M',values,'TolX',1e-3)

      use_lsqr = get_option(varargin,'lsqr',false);
      lsqr_tol = get_option(varargin,'lsqr_tol',1e-2);
      lsqr_iters = get_option(varargin,'lsqr_iters',50);
      if use_lsqr
          [w,flag,relres,iter,resvec] = lsqr(M',values,lsqr_tol,lsqr_iters);
          while flag > 0
             tolerance=tolerance*1.3; 
            [w,flag,relres,iter,resvec] = lsqr(M',values,tolerance,50);  
          end 
      else
          options = optimoptions('lsqlin','Algorithm','interior-point','Display','iter','TolCon',1e-10,'TolX',1e-14,'TolFun',1e-10);
          [n1,n2] = size(M');
          [w,resnorm,residual,exitflag,output,lambda_lsqlin] = ...
           lsqlin(M',values,-eye(n2,n2),zeros(n2,1),[],[],[],[],[],options);
      end
      
      ODF_stats = get_option(varargin,'ODFstats',false);
      if ODF_stats
         err = abs(M'*w - values);
         disp(['   Minimum weight: ',xnum2str(min(w))]);
         disp(['   Maximum weight:  : ',xnum2str(max(w))]);
         disp(['   Maximum error during interpolating ODF: ',xnum2str(max(err))]);
         disp(['   Mean error during interpolating ODF   : ',xnum2str(mean(err))]);
      end

      odf = odf + (1-m).*unimodalODF(S3G,psi,'weights',w./sum(w));
          
    end
    
    function odf = ambiguity1(varargin)

      cs = crystalSymmetry('222');

      orix = orientation.byAxisAngle(xvector,90*degree,cs);
      oriy = orientation.byAxisAngle(yvector,90*degree,cs);
      oriz = orientation.byAxisAngle(zvector,90*degree,cs);

      odf = unimodalODF([orix,oriy,oriz],varargin{:});
    end

    function odf = ambiguity2(varargin)
      cs = crystalSymmetry('222');
      ori = orientation.byAxisAngle(vector3d(1,1,1),[0,120,240]*degree,cs);
      odf = unimodalODF(ori,varargin{:});
    end
    
  end
     
end
