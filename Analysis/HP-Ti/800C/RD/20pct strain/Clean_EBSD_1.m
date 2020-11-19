%% Import the Data
% create an EBSD variable containing the data

if strcmp(opt.EBSD.interface,'mat')
    tmp=load(opt.EBSD.fname);
    ebsd=tmp.rebsd;
    ebsd=ebsd.gridify;
    ebsd.CSList=opt.CS;
elseif or(strcmp(opt.EBSD.interface,'osc'),strcmp(opt.EBSD.interface,'ang'))
    ebsd = loadEBSD(opt.EBSD.fname,opt.CS,'interface',opt.EBSD.interface,opt.EBSD.frame);
else
    error('unsupported file format for EBSD')
end

%% Clean Data
%color scheme of ipf maps
oM = ipfHSVKey(ebsd('indexed'));
% figure;plot(oM)

%Cleaned data
ebsd=ebsd(ebsd.prop.confidenceindex>opt.grain_recon.min_CI);

%Plot hex
color = oM.orientation2color(ebsd('indexed').orientations);
figure; plot(ebsd,color)
print('EBSD_Clean','-dtiffn','-r400') ;

%Plot grid
ebsdg=ebsd.gridify
color = oM.orientation2color(ebsdg('indexed').orientations);
figure; plot(ebsdg,color)
print('EBSDG_Clean','-dtiffn','-r400');

%Save cleaned EBSD
save('EBSDG_Clean.mat', 'ebsdg');
save('EBSD_Clean.mat', 'ebsd');


%% Crop small section 


ind = inpolygon(ebsd,[400,400,800,800]); % select indices by rectangle
ebsd_crop=ebsd(ind);
% ebsd_cropgrid=ebsd_crop.gridify;
oM = ipfHSVKey(ebsd('indexed'));

color = oM.orientation2color(ebsd_crop('indexed').orientations);
figure; plot(ebsd_crop,color)
print('EBSD_crop_Clean','-dtiffn','-r400');
% save('EBSDG_crop_Clean.mat', 'ebsd_cropgrid');
save('EBSD_crop_Clean.mat', 'ebsd_crop');