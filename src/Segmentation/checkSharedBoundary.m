function ShareBoundaryRatio = SharedBoundaryRatio(grains)
% check for points or grains to be inside a big grain
%
% Syntax
%   pctSharedBoundary = SharedBoundaryRatio(grains)
%
% Input
%  grains      - @grain2d

%
% Output
%  pctSharedBoundary - grains x grains matrix


% --- share boundary ratio computation---
nGrains=length(grains)
ShareBoundaryRatio = zeros(nGrains,nGrains);
mineral=grains.mineral; %For defining phase boundary in the case there is unindexed data
gBId=grains.boundary(mineral,mineral).grainId
gId = grains.id;
for i = 1:nGrains
   gBPairs=gBId(any(gBId==gId(i),2),:)
   [ugBPairs, ~, igBPairs] = unique(gBPairs, 'rows');
   ind=sub2ind([nGrains nGrains],ugBPairs(:,1),ugBPairs(:,2));
   ShareBoundaryRatio(ind) = accumarray(igBPairs, 1)./length(gBPairs);
end
