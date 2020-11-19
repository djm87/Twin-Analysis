%% Reconstruct and clean grains
%Calculate all grains
load('EBSDG_Clean.mat')
load('EBSD_Clean.mat')

[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',opt.grain_recon.seg_angle);

save('Grains_Clean.mat', 'grains');

figure;plot(grains,grains.meanOrientation)
print('grains','-dtiffn','-r400');

%Remove grains that are too small
grainsSelected = grains(grains.grainSize > opt.grain_recon.min_grainSz);
grainsSelectedx2 = grains(grains.grainSize > opt.grain_recon.min_grainSz*2);
grainsSelectedx3 = grains(grains.grainSize > opt.grain_recon.min_grainSz*3);

%reindex ebsd 
ebsdx3 = ebsd(grainsSelectedx3);
ebsdx2 = ebsd(grainsSelectedx2);
ebsdx1 = ebsd(grainsSelected);

% and perform grain reconstruction with the reduced EBSD data set
[grainsx1,ebsdx1.grainId,ebsdx1.mis2mean] = calcGrains(ebsdx1,'angle',opt.grain_recon.seg_angle);
[grainsx2,ebsdx2.grainId,ebsdx2.mis2mean] = calcGrains(ebsdx2,'angle',opt.grain_recon.seg_angle);
[grainsx3,ebsdx3.grainId,ebsdx3.mis2mean] = calcGrains(ebsdx3,'angle',opt.grain_recon.seg_angle);

figure;plot(grainsx1,grainsx1.meanOrientation)
print('grainsx1','-dtiffn','-r400');
figure;plot(grainsx2,grainsx2.meanOrientation)
print('grainsx2','-dtiffn','-r400');
figure;plot(grainsx3,grainsx3.meanOrientation)
print('grainsx3','-dtiffn','-r400');

%Save cleaned EBSD
grains=grainsx1;ebsd=ebsdx1;
save('Grainsx1_Clean.mat', 'grains');
save('EBSDx1_Clean.mat', 'ebsd');
grains=grainsx2;ebsd=ebsdx2;
save('Grainsx2_Clean.mat', 'grains');
save('EBSDx2_Clean.mat', 'ebsd');
grains=grainsx3;ebsd=ebsdx3;
save('Grainsx3_Clean.mat', 'grains');
save('EBSDx3_Clean.mat', 'ebsd');

%% Experiment reconstruction and clean grain
%Calculate all grains
% load('EBSDG_crop_Clean.mat')
Inputs;
load('EBSD_crop_Clean.mat')
F = meanFilter;
F = splineFilter;

% define the size of the window to be used for finding the median
% F.numNeighbours = 5; % this corresponds to a 7x7 window

[grains,ebsd_crop.grainId,ebsd_crop.mis2mean] = calcGrains(ebsd_crop,'angle',6*degree);
grainsSelected = grains(grains.grainSize > opt.grain_recon.min_grainSz*3);

%reindex ebsd 
ebsd_crop = ebsd_crop(grainsSelected);
ebsdS = smooth(ebsd_crop,F);
ebsdS = ebsdS('indexed');

% and perform grain reconstruction with the reduced EBSD data set
[grains,ebsdS.grainId,ebsdS.mis2mean] = calcGrains(ebsdS,'angle',6*degree);

grains=smooth(grains)

figure;plot(grains,grains.meanOrientation)

ebsd=ebsdS;
save('Grains_crop_Clean.mat', 'grains');
save('EBSD_crop_Clean.mat', 'ebsd');

