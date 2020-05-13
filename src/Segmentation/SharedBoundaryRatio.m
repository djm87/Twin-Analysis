function SBR = SharedBoundaryRatio(grains)
% finds the ratio of boundary shared by a grain with its neighbors
%
% Syntax
%   SharedBoundaryRatio = SharedBoundaryRatio(grains)
%
% Input
%  grains      - @grain2d

%
% Output
%  SharedBoundaryRatio - grains x grains matrix


% --- share boundary ratio computation---

% Initialize array
nGrains=length(grains);
SBR = zeros(nGrains,nGrains);

% Extract grainboundary grain ids
% mineral=grains.mineral; %For defining phase boundary in the case there is unindexed data
% gBId=grains.boundary(mineral,mineral).grainId;
gBId=grains.boundary.grainId;

% Get the grain Ids
gId = grains.id;

% Loop over grains computing the shared boundary ratios
for i = 1:nGrains
   % Find the unique gB pairs
   gBPairs=gBId(any(gBId==gId(i),2),:);
   % flip so grain looking at is in first column
   ind=gBPairs==gId(i);
   gBPairs(ind(:,2),:)=fliplr(gBPairs(ind(:,2),:));
   [ugBPairs, ~, igBPairs] = unique(gBPairs, 'rows');
   
   %Remove the boundary phase from the indexing
   if any(any(ugBPairs==0,2))
       ugBPairs(any(ugBPairs==0,2),:)=[];
       igBPairs(any(gBPairs==0,2))=[];
       igBPairs=igBPairs-1;
   end
   
   if ~isempty(ugBPairs)
       % Get linear index for ratio in SharedBoundaryRatio
       ind=sub2ind([nGrains nGrains],ugBPairs(:,1),ugBPairs(:,2));
    
       % Compute ratios and store
        
       SBR(ind) = accumarray(igBPairs, 1)./size(gBPairs,1);
       %    [accumarray(igBPairs, 1)./length(gBPairs),ugBPairs(:,2)]
   end
end
