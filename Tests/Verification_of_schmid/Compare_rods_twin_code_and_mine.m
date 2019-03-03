basedir='E:\Dropbox\Research\Source\Twin-Analysis\Rods Data\Zr-LN-IP-5\LN-IP-05-5\'

%% Read in the frags
frags.fname='frag.xlsx'
frags.sheet = 1;
frags.xlRangefragId = 'B2:B576';
frags.xlRangeGrainId = 'H2:H576';
frags.xlRangeType = 'J2:J576';
frags.xlRangeGen = 'K2:K576';
frags.xlRangeEffSchmid = 'N2:N576'; 

frags.fragId = xlsread(frags.fname,frags.sheet,frags.xlRangefragId)
frags.GrainId = xlsread(frags.fname,frags.sheet,frags.xlRangeGrainId)
frags.Type = xlsread(frags.fname,frags.sheet,frags.xlRangeType)
frags.Gen = xlsread(frags.fname,frags.sheet,frags.xlRangeGen)
frags.EffSchmid = xlsread(frags.fname,frags.sheet,frags.xlRangeEffSchmid)

%% Read in the grains
grains.fname='grains.xlsx'
grains.sheet = 1;
grains.xlRangeEuler = 'D2:F211';
grains.xlRangeId = 'B2:B211'; 

grains.euler = xlsread(grains.fname,grains.sheet,grains.xlRangeEuler)
grains.id = xlsread(grains.fname,grains.sheet,grains.xlRangeId)

%% Compute per grain schmid factors
CS = {'notIndexed',crystalSymmetry('622', [3.2 3.2 5.2],...
    'X||a', 'Y||b*', 'Z||c', 'mineral', 'Ti', 'color', 'light blue')};

CRSS=[100,133,125,125,111,111,111,111];

sS={slipSystem(Miller(1,0,-1,-1,CS{2},'uvtw'), Miller(1,0,-1,2,CS{2},'hkl'),CRSS(1)),... %T1
    slipSystem(Miller(-1,-1,2,6,CS{2},'uvtw'), Miller(1,1,-2,1,CS{2},'hkl'),CRSS(1)),... %T2
    slipSystem(Miller(1,1,-2,-3,CS{2},'uvtw'), Miller(1,1,-2,2,CS{2},'hkl'),CRSS(1))};    %C1

%apply symmetries for full slip systems
% for i=1:length(sS) 
% %     sS{i}=sS{i}.symmetrise('antipodal');
% end

% sigma = tensor([0 0 0; 0 0 0; 0 0 -1],'name','stress') %Sign of loading does more than just invert twin/parent flag
sigma = -stressTensor.uniaxial(vector3d.Z)

%% Find the computed schmid factors
computed=and(frags.EffSchmid~=0,~isnan(frags.EffSchmid))
EffSchmidRods=frags.EffSchmid(computed)
rawPhi1=grains.euler(:,1);
rawPHI=grains.euler(:,2);
rawPhi2=grains.euler(:,3);
ori=orientation('Euler',rawPhi1,rawPHI,rawPhi2,CS)
ori=orientation('Euler',3.64441,	0.0232034,	2.65164,CS)
rCS1=rotate(sigma,inv(ori));


[tauMax,m,n,tau1] = calcShearStress(rCS1,sS{1}.n,sS{1}.b,'symmetrise')
[tauMax,m,n,tau2] = calcShearStress(rCS1,sS{2}.n,sS{2}.b,'symmetrise')
[tauMax,m,n,tau3] = calcShearStress(rCS1,sS{3}.n,sS{3}.b,'symmetrise')
tau1=tau1';
tau2=tau2';
tau3=tau3';
tau1
% grainId=frags.GrainId(computed)

%Look at element 
i=1
tau2
%%
[SFMin1,active1] = max((SF1),[],2);
[SFMin2,active2] = max((SF2),[],2);
[SFMin3,active3] = max((SF3),[],2);
% SFMin1=SF1(active1);
% SFMin2=SF2(active2);
% SFMin3=SF3(active3);
SF=zeros(length(ori),1);
SF(frags.Type(computed)==0)=0;
SF(frags.Type(computed)==2)=SFMin1(frags.Type(computed)==2);
SF(frags.Type(computed)==3)=SFMin2(frags.Type(computed)==3);
SF(frags.Type(computed)==4)=SFMin3(frags.Type(computed)==4);

% SF=SFMin1
[EffSchmidRods,SF]

%% 

CS = crystalSymmetry('cubic',[3.523,3.523,3.523],'mineral','Nickel')
n = Miller(1,1,1,CS,'hkl')
d = Miller(0,-1,1,CS,'uvw')
r = normalize(vector3d(1,2,3))
tau = dot(d,r,'noSymmetry') * dot(n,r,'noSymmetry')
sS = slipSystem(d,n)
SF2 = sS.SchmidFactor(o\stressTensor.uniaxial(v));