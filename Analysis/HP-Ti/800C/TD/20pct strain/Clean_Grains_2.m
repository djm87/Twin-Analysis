%% Load the inputs
Inputs

%% Reconstruct and clean grains
%Calculate all grains
load('EBSDG_Clean.mat')
load('EBSD_Clean.mat')

[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',seg_angle);

save('Grains_Clean.mat', 'grains');

figure;plot(grains,grains.meanOrientation)
print('grains','-dtiffn','-r400');

%Remove grains that are too small
grainsSelected = grains(grains.grainSize > min_grainSz);
grainsSelectedx2 = grains(grains.grainSize > min_grainSz*2);
grainsSelectedx3 = grains(grains.grainSize > min_grainSz*3);

%reindex ebsd 
ebsdx3 = ebsd(grainsSelectedx3);
ebsdx2 = ebsd(grainsSelectedx2);
ebsdx1 = ebsd(grainsSelected);

% and perform grain reconstruction with the reduced EBSD data set
[grainsx1,ebsdx1.grainId,ebsdx1.mis2mean] = calcGrains(ebsdx1,'angle',seg_angle);
[grainsx2,ebsdx2.grainId,ebsdx2.mis2mean] = calcGrains(ebsdx2,'angle',seg_angle);
[grainsx3,ebsdx3.grainId,ebsdx3.mis2mean] = calcGrains(ebsdx3,'angle',seg_angle);

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