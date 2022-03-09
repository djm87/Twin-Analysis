%% Load and save EBSD locally 
% CS = {'notIndexed',crystalSymmetry('222', [2.840 5.870 4.940], 'mineral', 'Uranium', 'color', 'light blue')};
% fname='C:\Users\dansa\Documents\MATLAB\Twin-Analysis\Data\Rods Data\U-EWC5\2010-01-15\Scan1\scan1_cleaned.osc';
% ebsd = EBSD.load(fname,CS,'interface','osc',...
%             'convertEuler2SpatialReferenceFrame','setting 2');
% save('U-EWC5.mat','ebsd')
%% Reconstruct grains
load('U-EWC5.mat');

seg_angle=3*degree;
min_grain_sz=4;

[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'angle',seg_angle);

% remove 4 pixel grains and any grains with nan orientations
ebsd(grains(grains.grainSize<=min_grain_sz)) = [];
ebsd(grains(isnan(grains.meanOrientation))) = [];
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'angle',seg_angle);

setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outofPlane');

figure;plot(grains,grains.meanOrientation)
figure;plot(ebsd,ebsd.orientations)
% figure;plot(ebsd,ebsd.prop.imagequality)
%% get weight of fragments in initial ODF 
% tmpODF=load('ND_800C_initial_rodf.mat');
% initialODF=tmpODF.rodf;
% 
% opt.ncpu=6;
% poolobj = gcp('nocreate'); % If no pool, do not create new one.
% if isempty(poolobj)
%     poolobj=parpool(opt.ncpu)
%     poolobj.IdleTimeout=30;
% end
% 
% radius=10*degree
% volumeInitialODF=zeros(length(grains),1);
% ori=grains.meanOrientation;
% parfor i=1:length(grains)
% %     volumeInitialODF5(i) = volume(initialODF,grains.meanOrientation(i),radius5)
%     volumeInitialODF(i) = volume(initialODF,ori(i),radius);
% end
% 
% grains.prop.volInitialODF=abs(round(volumeInitialODF,6,'significant'));
% figure;plot(grains,grains.prop.volInitialODF)     


%% Save the ebsd and grains for use in segmentation
save('grains_recon.mat','grains')
save('ebsd_recon.mat','ebsd')