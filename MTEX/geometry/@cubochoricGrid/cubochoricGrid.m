classdef (InferiorClasses = {?rotation,?quaternion}) cubochoricGrid < orientation
%  
% Syntax
%  S3G = cubochoricGrid(ori,alphabeta,gamma)
%
% Input
%  points     - number of nodes
%  nodes      - @quaternion
%  resolution - double
%  CS, SS     - @symmetry groups
%
% Options
%  regular    - construct a regular grid
%  equispaced - construct a equispaced grid%
%  phi        - use phi
%  ZXZ, Bunge - Bunge (phi1 Phi phi2) convention
%  ZYZ, ABG   - Matthies (alpha beta gamma) convention
%  maxAngle   - only up to maximum rotational angle
%  center     - with respect to this given center

  
  properties
    
  end
  
  methods
    function S3G = cubochoricGrid(x,y,z,varargin)
      
      % call superclass method
      S3G = S3G@orientation(ori);

    
    end
  end
end
