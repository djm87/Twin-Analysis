% load some example data
mtexdata twins silent

% segment grains
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),...
    'angle',5*degree);

% remove two pixel grains
ebsd(grains(grains.grainSize<=2)) = [];
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),...
    'angle',5*degree,'removeQuadruplePoints');

%smooth grains
grains = grains.smooth(5);

%plot map
figure; plot(grains,grains.meanOrientation);hold on
plot(grains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-',...
  'displayName','grain boundary')
saveFigure('as_imported.png')

CS=grains.CS
k1=Miller(1,0,-1,2,twin{tnum}.CS,'hkl'); %What you specify here affects sign and schmid value
eta1=Miller(-1,0,1,1,twin{tnum}.CS,'uvtw'); %What you specify here affects sign and schmid value

%Define the Rtw such that Rtw transforms crystal to twin or
%phase to phase if not twin
mori{i}.Rtw(j)=orientation.map(mori{i}.u1,mori{i}.v1,...
    mori{i}.u2,mori{i}.v2); 

twinning = orientation.map(Miller(0,1,-1,-2,CS),Miller(0,-1,1,-2,CS),...
  Miller(2,-1,-1,0,CS),Miller(2,-1,-1,0,CS));
round(twinning.axis)
twinning.angle/degree
% extract all Magnesium Magnesium grain boundaries
gB = grains.boundary('Magnesium','Magnesium');

% and check which of them are twinning boundaries with threshold 5 degree
isTwinning = angle(gB.misorientation,twinning) < 5*degree;
twinBoundary = gB(isTwinning)

plot(twinBoundary,'linecolor','w','linewidth',5,'displayName',...
    'twin boundary');hold off
saveFigure('with_twin_boundaries.png')

%Merge twin boundaries
[mergedGrains,parentId] = merge(grains,twinBoundary);

% plot merged map
figure; plot(grains,grains.meanOrientation);hold on
plot(mergedGrains.boundary,'linecolor','k','linewidth',10,'linestyle','-',...
  'displayName','merged grains');
hold off
saveFigure('merged.png')



%% Plot texture 
%Plot c-axis distribution
figure
plotPDF(grains.meanOrientation,{Miller(0,0,0,1,CS)},'smooth')
saveFigure('texture_total.png')

%Select c-axis that center around the zvector and plot
cdir=Miller(0,0,0,1,CS)
grain_cdir=grains.meanOrientation.*cdir
lid = angle(grain_cdir, zvector)/degree < 45 | ...
    angle(grain_cdir,-zvector)/degree < 45

figure;plotPDF(grains(lid).meanOrientation,{Miller(0,0,0,1,CS)},'smooth')
saveFigure('texture_center.png')


%Plot twin (1) and non-twin (0)
figure;plot(grains,int8(lid))
hold on
plot(mergedGrains.boundary,'linecolor','k','linewidth',2.5,'linestyle','-')
hold off
saveFigure('texture_twin_map.png')

%% Schmid type of things 
 % Specimen stress state
 opt.sigma=stressTensor([0 0 0; 0 0 0; 0 0 -1])


 
 
 
 