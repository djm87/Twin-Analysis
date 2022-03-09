mtexdata twins silent
seg_angle=3*degree;
min_grain_sz=2;

[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'angle',seg_angle);

% remove two pixel grains
ebsd(grains(grains.grainSize<=min_grain_sz)) = [];
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'angle',seg_angle);

% smooth the grains
grains = grains.smooth(3);

save('grains_recon.mat','grains')
save('ebsd_recon.mat','ebsd')