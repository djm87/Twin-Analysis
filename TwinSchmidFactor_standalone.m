%% General inputs
% Note that to simplify things, the schmid factor output will always be >0
% for slip and this code does not distinguish if the slip is in the
% negative direction.  For twinning, schmid factor output will always be
% -90 to +90, as it is directional. note that for magnesium a plane and the
% normal will not necessarily have the same index values (unlike cubic
% systems), so use the ,'hkl' suffix for plane normals and the 'uvw' suffix
% for directions. %

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');

%define system crystallography
CS = {... 
  'notIndexed',...
  crystalSymmetry('622', [3.2 3.2 5.2], 'X||a', 'Y||b*', 'Z||c', 'mineral', 'Magnesium', 'color', 'light blue')};
% CS{2} = crystalSymmetry('622', [3.2 3.2 5.2], 'X||a', 'Y||b*', 'Z||c', 'mineral', 'Magnesium', 'color', 'light blue');



% path to files
pname = 'E:\Dropbox\Research\Tutorials\MTEX Twin Area extraction\WE43\0.08strain';

% which files to be imported
fname = [pname '\0.08strain.ang'];

% create an EBSD variable containing the data
ebsd = loadEBSD(fname,CS,'interface','ang',...
  'convertEuler2SpatialReferenceFrame');
oM = ipdfHSVOrientationMapping(ebsd('indexed'));
color = oM.orientation2color(ebsd('indexed').orientations);
% figure;plot(oM); %IPF triangle
figure;plot(ebsd('indexed'),color); %EBSD IPF map

%define tensile axis or for an arbitrary vector use r = vector3d(1,-1,0);
r=xvector; %not x

%% define slip systems and twinnings
% define  CRSS of the basal, prismatic, pyramidal, extension twin so that
% the schmid factor can be scaled relatively between the systems based on
% Figure 11, Chapuis, A. & Driver, J. H. Temperature dependency of slip and
% twinning in plane strain compressed magnesium single crystals Acta
% Mater., 2011, 59, 1986-1994 in MPa.  BasalFactor is 4.  Either input
% directly
%Set for WE43
sSBasal = slipSystem(Miller(2,-1,-1,0,CS{2},'uvtw'), Miller(0,0,0,1,CS{2},'hkl'),20/30);
sSBasal = sSBasal.symmetrise%('antipodal')

sSPrismatic = slipSystem(Miller(-1,2,-1,0,CS{2},'uvtw'), Miller(1,0,-1,1,CS{2},'hkl'),70/30);
sSPrismatic = sSPrismatic.symmetrise%('antipodal')

% type II pyramidal
sSPyramidal = slipSystem(Miller(2,-1,-1,-3,CS{2},'uvtw'), Miller(2,-1,-1,2,CS{2},'hkl'),215/30);
sSPyramidal = sSPyramidal.symmetrise%('antipodal')

%{10-12} extension twinning
sSExtTwin1 = slipSystem(Miller(-1,0,1,1,CS{2},'uvtw'), Miller(1,0,-1,2,CS{2},'hkl'),111/30);
sSExtTwin1 = sSExtTwin1.symmetrise%('antipodal')

%{11-21} extension twinning
sSExtTwin2 = slipSystem(Miller(1,1,-2,6,CS{2},'uvtw'), Miller(-1,-1,2,1,CS{2},'hkl'),113/30);
sSExtTwin2 = sSExtTwin2.symmetrise%('antipodal');

%{10-11} contraction twinning
sSConTwin1 = slipSystem(Miller(2,-1,-1,-3,CS{2},'uvtw'), Miller(2,-1,-1,2,CS{2},'hkl'),195/30);
sSConTwin1 = sSConTwin1.symmetrise%('antipodal');

%% Compute grains to get mean orientation
% segmentation angle typically 10 to 15 degrees that seperates to grains
seg_angle = 5;

% minimum indexed points per grain between 5 and 10
min_points = 10;

%Restrict to indexed points
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),...
    'angle',seg_angle*degree);


% tensor direction in crystal coordinates
% rCS = ebsd('Magnesium').orientations \ r;
rCS = grains.meanOrientation \ r;
%% Basal slip calculations

sFBasal = abs(sSBasal.SchmidFactor(rCS));

% figure;plot(ebsd('Magnesium'),max(sFBasal,[],2))

%% Prismatic slip calculations

sFPrismatic = abs(sSPrismatic.SchmidFactor(rCS));
% figure;plot(ebsd('Magnesium'),max(sFPrismatic,[],2))

%% Pyramidal slip calculations

sFPyramidal = abs(sSPyramidal.SchmidFactor(rCS));
% figure;plot(ebsd('Magnesium'),max(sFPyramidal,[],2))

%% Ext Twin 1 calculations

sFExtTwin1 = abs(sSExtTwin1.SchmidFactor(rCS));
% figure;plot(ebsd('Magnesium'),max(sFExtTwin1,[],2))

%% Ext Twin 2 calculations

sFExtTwin2 = abs(sSExtTwin2.SchmidFactor(rCS));
% figure;plot(ebsd('Magnesium'),max(sFExtTwin2,[],2))

%% Contraction Twin 1 calculations

sFConTwin1 = abs(sSConTwin1.SchmidFactor(rCS));
% figure;plot(ebsd('Magnesium'),max(sFConTwin1,[],2))

%% compare SchmidFactors scaled by CRSS

% combine all slip systems
sS = [sSBasal;sSPrismatic;sSPyramidal;sSExtTwin1;sSExtTwin2;sSConTwin1];

% compute the relative Schmid factors scaled by CRSS
sFRelative = sS.SchmidFactor(rCS,'relative');

%% 

% compute the maximum relative Schmid factors
[sFmax,sFid] = max(abs(sFRelative),[],2);

sFRelativeMaxMag=zeros(length(sFmax),1);
for i=1:length(sFmax)
    sFRelativeMaxMag(i)=sFRelative(i,sFid(i));
end

% plot the relatibe 
figure(1); plot(grains,double(sFRelativeMaxMag>0))

%to map the active slip or twinning system
% figure(2);
% 
% plot(ebsd(sFid > 0 & sFid <= 6),'facecolor',[0.5 0.5 1],'DisplayName','basal')
% hold on
% plot(ebsd(sFid > 6 & sFid <= 12),'facecolor',[0.5 1 0.5],'DisplayName','prismatic')
% plot(ebsd(sFid > 12 & sFid <= 24),'facecolor',[1 0.5 0.5],'DisplayName','pyramidal')
% plot(ebsd(sFid > 24 & sFid <= 36),'facecolor','yellow','DisplayName','ext twin1')
% plot(ebsd(sFid > 36 & sFid <= 24),'facecolor','yellow','DisplayName','ext twin2')
% plot(ebsd(sFid > 24 & sFid <= 30),'facecolor','yellow','DisplayName','con twin')

% hold off