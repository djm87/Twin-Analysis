%% CRYSTALLOGRAPHIC VARIABLE LIST
%This code analyses EBSD data and outputs information for all
%mergedGrains and grains.  This version of the code assumes the largest 
%grain composing the merged grain is the parent.  For manual selection of 
%the parent see Parent_Twin_Interact instead.
%Customization options include material, twinning type,
%and direction of applied force (r)

%input; EBSD data (ebsd variable).

%output; two arrays: TwAr which includes information on the mergedGrains
%such as area, area of twins (limited to those over the minSize), count of
%twins, broken down by variant.  Other output is GrAr which includes one line
%per grain with information on the size, aspect ratio, and labels it 'P'
%for parent, 'T' for twin, 'A' for secondary twin and 'O' for other.  All
%SF (Schmid factors) in these are those for twinning only.

%Customization; This code assumes you are dealing with magnesium extension twinning.
%To use the code for another material, replace 'Magnesium' as appropriate
%below, and input the correct twinning parameters below.  All standard
%customization options are in the current section

%R7 - updated for consistency with changes made to reload code.  Incoprporates
%additional check for secondary twins labeled A and summarises counts and
%areas. 

%R8 fixed math error in calculation of major/minor axis.  Improved comments

%to run this code, EBSD data must be loaded in the workspace as a variable
%import the symmetry
cs=ebsd('Magnesium').CS;
%for schmid calculations, enter the vector3d along which the force is
%applied (default is x)
r=xvector;

%define the minimum pixel count below which grain is not considered the
%parent
minSize=8;


% n1 – the twinning direction or the direction of shear lying in K1; for Mg extension
% twinning this is family of planes described by Miller(1,0,-1,1,cs,'uvw')
% note that 'uvw' shows the indices are of a vector. refered to in the code
%as En1 (for extension twin n direction)
En1 = [Miller(1,0,-1,1,cs,'uvw');  Miller(-1,0,1,1,cs,'uvw');Miller(-1,1,0,1,cs,'uvw');Miller(1,-1,0,1,cs,'uvw');Miller(0,-1,1,1,cs,'uvw');Miller(0,1,-1,1,cs,'uvw')];

%(1) K1 – the twinning or composition plane that is the invariant
% (unrotated and undistorted) plane of the simple shear; for Mg extension
% twinning this is family of planes described by Miller(-1,0,1,2,cs,'hkl').
% note that 'hkl' shows the indices are of a plane normal which is critical
% in HCP as the directions and the planes are not necessarily
% perpendicular.  These are refered to in the code as Ek1 (for extension
% twin K1 plane)
%Full list of all  normals to the slip planes
Ek1=[ Miller(-1,0,1,2,cs,'hkl');  Miller(1,0,-1,2,cs,'hkl'); Miller(1,-1,0,2,cs,'hkl');Miller(-1,1,0,2,cs,'hkl'); Miller(0,1,-1,2,cs,'hkl');Miller(0,-1,1,2,cs,'hkl')];

% Extension twin Axis around which the twin rotates
Eta=[Miller(-1,2,-1,0,cs,'uvw');Miller(1,-2,1,0,cs,'uvw');Miller(-1,-1,2,0,cs,'uvw');Miller(1,1,-2,0,cs,'uvw');Miller(2,-1,-1,0,cs,'uvw');Miller(-2,1,1,0,cs,'uvw')];

%define an extension twin (86.29 degrees reorientation) by axis and angle of rotation
%this is used for detecting twin boundaries
h= Miller(1,1,-2,0, cs);
%for extension twin
twinning = orientation('axis',h,'angle',86.3*degree,cs,cs);


%% Threshold into grains, and merge twins into parents

%threshold into grains
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'angle',10*degree);

%Extract the GB information to a separate variable
gB = grains.boundary;

%To not include the grain-border in calculations, create a new variable for only the mg-mg segments
gB_MgMg = gB('Magnesium','Magnesium');

%Next we check for each boundary segment whether it is a twinning boundary, %i.e., whether boundary misorientation is close to the twinning. 
% restrict to twinnings with threshold 5 degree
isTwinning = angle(gB_MgMg.misorientation,twinning) < 5*degree;
twinBoundary = gB_MgMg(isTwinning);

[mergedGrains,parentId] = merge(grains,twinBoundary);

%% Set up a matrix with one line per mergedGrain

%name all the columns
columntitle = {'merged ID','merged parent area','M.P. major axis','M.P. aspect ratio','M.P. percent internal boundary length',...
'parent x','parent y','Orientation of parent at xy','GOS of main parent','Count Parent grains','Total Parent Area',...
'Count Twin grains','Total Twinned Area','Total Other Area','Percent of grain that twinned',...
'Schmid for Var1','Schmid for Var2','Schmid for Var3','Schmid for Var4','Schmid for Var5','Schmid for Var6',...
'Area Twinned for Var1','Area Twinned for Var2','Area Twinned for Var3','Area Twinned for Var4','Area Twinned for Var5','Area Twinned for Var6'...
'Rank for Var1','Rank for Var2','Rank for Var3','Rank for Var4','Rank for Var5','Rank for Var6',...
'Twin Count for Var1','Twin Count for Var2','Twin Count for Var3','Twin Count for Var4','Twin Count for Var5','Twin Count for Var6'...
'M.P. total boundary length','M.P. border boundary length','Parent grain','total secondary twins','total Area secondary twins'};
%This array is information on the mergedGrains, one row per entry.
%in order this contains the ID#, area, the longest length, the aspect
%ratio, the percent of boundary that does not contact the ebsd scan borders
%the x and y coordinates of the point used as the parent orientation, the
%euler angles (bunge convention) of that orientation, the grain orientation
%spread of the parent grain at x,y, the count of grains with orientaiton 
%within 10 degrees of parent that are larger than the min size, the area of
%these grains, the number of grains within 10 degrees of the possible twin 
%variant orientations that are larger than the min size, the area of
%these grains, the area of all other grains larger than the min size, the
%percent of the grain twinned (100*twin area /parent plus the twin area),
%schmid factor for all six twin variants, the area twinned for each
%variant, the rank of all six variants (highest SF=1, lowest =6), the count
%of twins for each variant type, total boundary length, and total boundary length
%along the border, the ID of the parent grain, and finally infromation
%about the secondary twins (provided > minSize)

%create an array with one line for each of the merged grains, and add the
%appropriate column headings
TwAr = cell(length(mergedGrains), length(columntitle)); 
TwAr = [columntitle;TwAr];
%% Set up a matrix with one line per grain

columntitle = {'grain ID','merged ID','grain type','grain area','major axis','minor axis',...
'percent grain boundary internal','deg dev from ideal','GOS','Twin Rank','Active Twin SF'};

GrowCount=length(parentId);
%must be total columns of GrAr-2 to allow for m and parentID
GrAr=cell(GrowCount,length(columntitle)-2);

%the row number is the grain ID so number sequentially
m = reshape(1:GrowCount, [1 GrowCount])';

%concatenate horizontally and add in the parentId
GrAr = [num2cell(m), num2cell(parentId), GrAr];

%label the columns by vertically concatenate with the titles 
GrAr = [columntitle; GrAr];

%% Loop to query parent xy, and retrieve info from grains

% determine the pixel size 
dx = max(ebsd.unitCell(:,1))-min(ebsd.unitCell(:,1));
dy = max(ebsd.unitCell(:,2))-min(ebsd.unitCell(:,2));
pixSize=abs(dx*dy);

%loop through all the merged parent data, once per entry in mergedGrains
for TwinNum=1:length(mergedGrains)
%output the data from the merged grain

%store the merged grain ID for reference
TwAr{(TwinNum+1),1}=mergedGrains(TwinNum).id;
%store the merged area NOTE if you use grain_selected.area this
%incorporates the GB curvature.  Better to use the # of pixels x pixel
%spacing.
area=full(mergedGrains(TwinNum).grainSize)*pixSize;
TwAr{(TwinNum+1),2}=area;

%store the grain major and minor axis.  Since we can easily get the area
%and the aspect ratio, we calculate based on area=pi*minor axis *major axis
%and minor axis * aspect ratio = major axis, so a(major axis/2)=
%area*aspect ratio divided by pi
majAxis=2*sqrt(area*(mergedGrains(TwinNum).aspectRatio)/pi);
%store the major axis length of the parent
TwAr{(TwinNum+1),3}=majAxis;
%store the aspect ratio
TwAr{(TwinNum+1),4}=mergedGrains(TwinNum).aspectRatio;

%store the total boundary length
boundZero=sum(any(mergedGrains(TwinNum).boundary.grainId==0,2));
boundTotal=mergedGrains(TwinNum).boundary.length;
boundPct=100*(boundTotal-boundZero)/boundTotal;
TwAr{(TwinNum+1),5}=boundPct;
TwAr{(TwinNum+1),40}=boundTotal;
TwAr{(TwinNum+1),41}=boundZero;

%Now, determine how many child grains went into this parent
%first put list of rows (which equal id # from grains variable) into an array
[a]=find(parentId==mergedGrains(TwinNum).id);

% the length of this tells how many entries from grains went into the merged grain
%entry (how many times this loop must run)
gLoop=length(a);

%Determine the parent based on which grain referenced by a is the largest
%first, get the area of all grains in a
sizea=full(grains(a).grainSize)*pixSize;
%determine the largest row
[M,I]=max(sizea);
%get the index in corresponding row of a
index = a(I);

%this grain is the parent.
TwAr{(TwinNum+1),42}=index;

%determine the parent orientation by selecting a point of the grain with
 %lowest angle to the meanOrientation, 
%first get id# for all EBSD orientations within that default parent grain
indexAll=find(ebsd.grainId==grains(index).id);
%now get the actual orientations included within default parent grain
oriGindex=ebsd(indexAll).orientations;
%now get the angle between the meanorientation of the parent and all other
%orientations within the parent grain
angG2M=angle(oriGindex, grains(index).meanOrientation)/degree;
%find the ebsd orientation that is closest to the meanOrientation
check1=abs(angG2M);
[angM1,angI1]=min(check1);
%get the x and y of this data
x=ebsd(indexAll(angI1)).x;
y=ebsd(indexAll(angI1)).y;

%Save these values of parent to the TwAr
TwAr{(TwinNum+1),6}=x;
TwAr{(TwinNum+1),7}=y;
oriP=ebsd(x,y).orientations;
[alpha,beta,gamma] = Euler(oriP,'bunge');
TwAr{(TwinNum+1),8} = strcat(num2str(rad2deg(alpha)),',',num2str(rad2deg(beta)),',',num2str(rad2deg(gamma)));

%record the GOS of the parent grain
TwAr{(TwinNum+1),9}=grains(index).GOS./degree;

%prealocate SF
SF=zeros(1,6);
%calculate the parent schmid factor 
for j=1:6
%calculate the angle of the twin plane normal in the parent
ExtTwinNormal= oriP*Ek1(j);
%angle between force and normal
theta=cos(angle(r, ExtTwinNormal));
%calculate the orientation of the slip directions in the parent
ExtTwinDir=oriP*En1(j);
%note that the 'antipodal' command is not used because twinning is
%directional.  this will give negative values of cos (range 1 to -1) when angle >90 degrees
%calculate the angle between the slip and the force
lmda=cos(angle(r,ExtTwinDir));
SF(j)=theta.*lmda;
% Now write schmid factor to array
colWrite=j+15;
TwAr{(TwinNum+1),colWrite}=SF(j);
end
%end of schmid factor calculation loop

%determine the rank of all schmid factors.  sortSF holds the SF in
%ascending order, and inv gives the rank in ascending order (highest schmid
%=6)
[sortSF, extra, invrank] = unique(SF);
%switch rank to descending order, ie highest schmid =1
rank=7-invrank;

%calculate the 6 possible orientations generated by this parent and record
%schmid rank too
for j=1:6
% Write schmid rank to array
colWrite=j+27;
TwAr{(TwinNum+1),colWrite}=rank(j);

%calculate the orientations of the six twin variants
%find the vector3d which is the direction of the twin axis in the parent
vari(j)=oriP*Eta(j);

%set a rotation around that axis
rot(j)=rotation('axis',vari(j),'angle',86.3*degree);

% rotate the c axis around the twin axis, and save orientation
oriV(j)=rot(j)*oriP;
end

%reset the running totals for number of twins and count
%also count of parent segments and running total of areas which don't match anything (ie errors)
totalAreaT=0;
totalCountT=0;
totalAreaP=0;
totalCountP=0;
totalAreaO=0;
totalCountA=0;
totalAreaA=0;

%these two are arrays of six rows (one per variant)
varCountT=zeros(1,6);
varAreaT=zeros(1,6);

%start loop here through all grains which went into merged grain 
for i=1:gLoop

%get the id of the grain to analyse in this loop from array a
index=a(i);

%and start storing data in GrAr, one line per grain
%and get the number of pixels
N=full(grains(index).grainSize);
%and save it as the area
GrAr{(index+1),4}=N*pixSize;

%store the grain major and minor axis.  Since we can easily get the area
%and the aspect ratio, we calculate based on area=pi*minor axis *major axis
%and minor axis * aspect ratio = major axis, so a(major axis/2)=
%area*aspect ratio divided by pi
majAxis=2*sqrt((N*pixSize)*(grains(index).aspectRatio)/pi);
GrAr{(index+1),5}=majAxis;
GrAr{(index+1),6}=majAxis/grains(index).aspectRatio;

%store the boundary info in the grain array
boundZero=sum(any(grains(index).boundary.grainId==0,2));
boundTotal=grains(index).boundary.length;
boundPct=100*(boundTotal-boundZero)/boundTotal;
GrAr{(index+1),7}=boundPct;

%record the orientation spread
GrAr{(index+1),9}=grains(index).GOS./degree;

%determine if the grain is part of the parent, a twin, or something else
ang=angle(oriP, grains(index).meanOrientation)/degree;
%check if it's a parent
if ang<9.999
%if it is, update GrAr with grain type
GrAr{(index+1),3}='P';
%and record deviation from ideal parent
GrAr{(index+1),8}=ang;

%if the is large enough, add to counters for parent # and %area
if N>minSize
totalCountP=totalCountP+1;
totalAreaP=totalAreaP+N*pixSize;
end
%and exit this loop of grain checking
continue
end

%if the grain is not a parent, we must check for twin boundary
%Extract the GB information to a separate variable
gB = grains(index).boundary;

%To not include the grain-border in calculations, create a new variable for 
%only the mg-mg segments
gB_MgMg = gB('Magnesium','Magnesium');

%Next we check for each boundary segment whether it is a twinning boundary, 
%i.e., whether boundary misorientation is close to the twinning. 
% restrict to twinnings with threshold 3 degree
isTwinning = angle(gB_MgMg.misorientation,twinning) < 3*degree;
%however, this is not a sufficient condition on its own because both sides
%of the boundary will qualify.  The twin must also be within a certain
%range of the allowable twin variants (this prevents 'other' type grains
%from qualifying erroniously)

%compare the 6 variant orientations to the grain in question, 
for j=1:6
angV(j)=angle(oriV(j), grains(index).meanOrientation)/degree;
end
%is the closest variant within 10 degrees?  be careful of negatives
check=abs(angV);
[angM,angI]=min(check);
%if the closest variant is within ~10 degrees of the grain in question,
%this is a twin.  
%only if thre exist a minimum of 3 boundary segments that qualify as twinning and
 %the orientation within 10 degrees of parent will the grain qualify as
 %twin
if sum(isTwinning)>3 && angM<9.999
%if it is, update GrAr
GrAr{(index+1),3}='T';

%and record deviation from ideal twin orientation
GrAr{(index+1),8}=angM;

%and use the index of the lowest angle to retrieve the corresponding schmid
%and rank
GrAr{(index+1),10}=rank(angI);
GrAr{(index+1),11}=SF(angI);

%if the grain has more than minSize, add to counters for twin # and
%area
if N>minSize
%add the area to the running total for twin area
totalAreaT=totalAreaT+N*pixSize;
%increment the twin count.
totalCountT=totalCountT+1;

%if it is a twin, the closest variant is in angI,
%increment the count for that variant

varCountT(angI)=varCountT(angI)+1;
%increment the area for that variant
varAreaT(angI)=varAreaT(angI)+N*pixSize;
%and exit this loop
end
continue
end

%if you're still in the loop, then it's something odd.  Increment area of
%other grains and finish loop if it's big enough to count
GrAr{(index+1),3}='O';
GrAr{(index+1),8}=ang;

if N>minSize
%if the grain is big enough, count the area
totalAreaO=totalAreaO+N*pixSize;
end

end
%end of loop through all grains in a given merged grain

%Check for secondary twins if both other grains and twins (of significant size) exist in
%the mergedGrain
if totalAreaO>0 && totalAreaT>0
%get a list of all applicable other grains and loop through them.
%list all other grains
typea=false(length(a),1);
for i=1:length(a)
row=a(i)+1;
if GrAr{row,3}=='O'
typea(i)=1;
else
typea(i)=0;
end
end
%eliminate all non-others from typea
outO=a(typea);

%for each 'other' grain determine if there is an adjacent twin with which it has a
%misorientation relationship by first checking the boundary orientation and
%relationship.
for i=1:length(outO)
%put the index of the applicable 'other grain' in variable j
indexA=outO(i);
bound_other=grains(indexA).boundary;
%To not include the grain-border in calculations, create a new variable for only the mg-mg segments
gB_MgMg = bound_other('Magnesium','Magnesium');
%Next we check for each boundary segment whether it is a twinning boundary, %i.e., whether boundary misorientation is close to the twinning. 
% restrict to twinnings with threshold 3 degree 

%creates a logical of which points are twins
isTwinning = angle(gB_MgMg.misorientation,twinning) < 3*degree;
%reduces the data set of boundaries to only include twins
twinBoundary = gB_MgMg(isTwinning);
%get the ID of all grains involved
CheckGrains=unique(twinBoundary.grainId);

%remove the 'O' grain from this list so we only have the grains on the
%opposite (outside) of the boundary.
CheckGrains(all(CheckGrains==indexA,2),:)=[];

%for each grain across a twin boundary, check if they are first within the grain and then
%classified as 'T'.  must be this order or may hit uncategorised grains on
%first run through of the original code.
for k=1:length(CheckGrains)
%first check if the grain is in a (which contains the list of all grains
%for this merged grain)
indexT=CheckGrains(k);
if any(a==indexT)
%now check if that grain is T
if GrAr{indexT+1,3}=='T'
%if it is, we've found a secondary twin.  Change the designation of the OTHER GRAIN to A
GrAr{indexA+1,3}='A';
%load both the primay twin (T) and secondary twin (A) orientations into variables.
oriT=grains(indexT).meanOrientation;
oriA=grains(indexA).meanOrientation;
%calculate the possible secondary twin variants 
%determine which one was created
%output the deviation from this orientation to GrAr column 8.

%prealocate SF and verify empty
SF=zeros(1,6);
angV=zeros(1,6);
clear vari
clear rot
clear oriV

%calculate the schmid factors of the secondary twin and the variant
%orientations
for j=1:6
%calculate the angle of the twin plane normal in the parent
ExtTwinNormal= oriT*Ek1(j);
%angle between force and normal
theta=cos(angle(r, ExtTwinNormal));
%calculate the orientation of the slip directions in the parent
ExtTwinDir=oriT*En1(j);
%note that the 'antipodal' command is not used because twinning is
%directional.  this will give negative values of cos (range 1 to -1) when angle >90 degrees
%calculate the angle between the slip and the force
lmda=cos(angle(r,ExtTwinDir));
SF(j)=theta.*lmda;

%calculate the orientations of the six twin variants
%find the vector3d which is the direction of the twin axis in the parent
vari(j)=oriT*Eta(j);

%set a rotation around that axis
rot(j)=rotation('axis',vari(j),'angle',86.3*degree);

% rotate the c axis around the twin axis, and save orientation of the
% variant
oriV(j)=rot(j)*oriT;
%save the deviation of the varinat from the other grain
angV(j)=angle(oriV(j), oriA)/degree;
end
%end of schmid factor calculation loop

check=abs(angV);
[angM,angI]=min(check);
%this is the deviation from an ideal secondary twin
GrAr{(indexA+1),8}=angM;

%determine the rank the active variant.  sortSF holds the SF in
%ascending order, and inv gives the rank in ascending order (highest schmid
%=6, opposite what we want)
[sortSF, extra, invrank] = unique(SF);
%switch rank to descending order, ie highest schmid =1
rank=7-invrank;
%and use the index of the lowest angle to retrieve the corresponding schmid
%and rank
GrAr{(indexA+1),10}=rank(angI);
GrAr{(indexA+1),11}=SF(angI);

%if the secondary twin is large enough add it to the counter and size running totals
N=full(grains(indexA).grainSize);
if N>minSize
totalCountA=totalCountA+1;
%calc area of other grain, and add to running total
totalAreaA=totalAreaA+N*pixSize;
end
%we can now proceed to the next 'other' grain on the list, no more twins
%need to be checked
break
else
continue
end
end
end
end
end
%Check for secondary twinning is now complete.

%write the running totals of perimiter and area to the
%matrix  We must remember to subtract the grains which were changed to
%secondary twins from the 'other' count
totalAreaO=totalAreaO-totalAreaA;
TwAr{(TwinNum+1),10}=totalCountP;
TwAr{(TwinNum+1),11}=totalAreaP;
TwAr{(TwinNum+1),12}=totalCountT;
TwAr{(TwinNum+1),13}=totalAreaT;
TwAr{(TwinNum+1),14}=totalAreaO;
TwAr{(TwinNum+1),43}=totalCountA;
TwAr{(TwinNum+1),44}=totalAreaA;
%only write other twin info if there are twins
if totalCountT>0
%calculate the percent of grain that twinned
TwAr{(TwinNum+1),15}=100*totalAreaT/(totalAreaP+totalAreaT);
%write twin areas by variant
TwAr{(TwinNum+1),22}=varAreaT(1);
TwAr{(TwinNum+1),23}=varAreaT(2);
TwAr{(TwinNum+1),24}=varAreaT(3);
TwAr{(TwinNum+1),25}=varAreaT(4);
TwAr{(TwinNum+1),26}=varAreaT(5);
TwAr{(TwinNum+1),27}=varAreaT(6);
%write twin counts by variant
TwAr{(TwinNum+1),34}=varCountT(1);
TwAr{(TwinNum+1),35}=varCountT(2);
TwAr{(TwinNum+1),36}=varCountT(3);
TwAr{(TwinNum+1),37}=varCountT(4);
TwAr{(TwinNum+1),38}=varCountT(5);
TwAr{(TwinNum+1),39}=varCountT(6);
end

end
%end loop through all merged grains
