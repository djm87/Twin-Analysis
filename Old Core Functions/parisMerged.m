function p = parisMerged(grains,CS)
% Percentile Average Relative Indented Surface
%
% the paris (Percentile Average Relative Indented Surface) is shape
% parameter that measures the convexity of a grain
%
% Syntax
%   p = paris(grains)
%
% Input
%  grains - @grain2d
%
% Output
%  p - double
%

p = zeros(size(grains));

% store this in local variables for speed reasons
X = grains.V(:,1);
Y = grains.V(:,2);

% poly = grains.poly;
% remove inclusions
% incl = grains.inclusionId;
% for i = find(incl>0).'
%   poly{i} = poly{i}(1:end-incl(i));
% end

% compute convex hull perimeters
for id = 1:length(grains)
  gB=grains(id).boundary;
  gB_mm=gB(CS.mineral,CS.mineral);
  % for small grains there is no difference
  if length(gB_mm) <= 7, continue; end
  
  % extract coordinates
  xGrain = gB_mm.x;%X(poly{id});
  yGrain = gB_mm.y;%Y(poly{id});
  
  % compute perimeter
  perimeterGrain = sum(sqrt(...
    (xGrain(1:end-1) - xGrain(2:end)).^2 + ...
    (yGrain(1:end-1) - yGrain(2:end)).^2));
   
  % compute convex hull
  ixy = convhull(xGrain,yGrain);
  
  % compute perimenter
  perimeterHull = sum(sqrt(...
    (xGrain(ixy(1:end-1)) - xGrain(ixy(2:end))).^2 + ...
    (yGrain(ixy(1:end-1)) - yGrain(ixy(2:end))).^2));
  
  % paris is the relative difference between convex hull perimenter and true
  % perimeter
  p(id) = 200*(perimeterGrain - perimeterHull)./perimeterHull;
  
end

