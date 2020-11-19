%% Import the Data
% create an EBSD variable containing the data
if strcmp(EBSDinterface,'mat')
    tmp=load(EBSDfname);
    ebsd=tmp.rebsd;
    ebsd=ebsd.gridify;
    ebsd.CSList=CS;
elseif or(strcmp(EBSDinterface,'osc'),strcmp(EBSDinterface,'ang'))
    ebsd = loadEBSD(EBSDfname,CS,'interface',EBSDinterface,EBSDframe);
else
    error('unsupported file format for EBSD')
end

%% Clean Data
%color scheme of ipf maps
oM = ipfHSVKey(ebsd('indexed'));
% figure;plot(oM)

%Cleaned data
ebsd=ebsd(ebsd.prop.confidenceindex>min_CI);

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
