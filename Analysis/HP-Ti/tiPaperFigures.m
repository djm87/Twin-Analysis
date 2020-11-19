%% Options
%clc
clear
%close all;

%% Parameters
start = 1;  %begin loop
fin = 12;    %end loop
ebsd_rot_angle = 0;

min_CI = 0.1;
min_grainSz = 10;
seg_angle = 5*degree;

%% Specify files to be imported, Crystal and Specimen Symmetries, and reference frames
% which files to be imported
baseDir=pwd;
temps={'625C','800C'}
direcitons={'ND','TD','RD'}
strain_lvls={'5pct strain','20pct strain'}
ebsdName='EBSDx3_Clean.mat'
grainName='Grainsx3_Clean.mat'
ebsdPath={}; fnameOut={}; grainPath={};
cnt=1;
for i=1:2
    for j=1:3
        for k=1:2
            ebsdPath{cnt}=fullfile(baseDir,temps{i},direcitons{j},strain_lvls{k},ebsdName)
            grainPath{cnt}=fullfile(baseDir,temps{i},direcitons{j},strain_lvls{k},grainName)
            fnameOut{cnt}=strrep(strcat(temps{i},'-',strain_lvls{k},'-',direcitons{j}),' ','-')
            cnt=cnt+1;
        end
    end
end

% crystal symmetry
CS = crystalSymmetry('622', [2.95 2.95 4.6855], 'X||a', 'Y||b*', 'Z||c',...
    'mineral', 'Titanium (Alpha)', 'color', 'light blue');

% EBSD frame
setMTEXpref('xAxisDirection','south');
setMTEXpref('zAxisDirection','outofPlane');

%Create file to write grain data to
fileID = fopen([baseDir,'\Figures\','Grain Size Data.txt'],'w');

%% Initiate loop
for num = start:fin   
    %% Load and rotate EBSD scan
    fprintf('Loading EBSD data for %s...\n',fnameOut{num});  %Display name of current EBSD scan

    load(ebsdPath{num});
    rebsd=ebsd;
    load(grainPath{num});

    grains =smooth(grains,4); % smooth the grains a bit
    rebsd_lowCI = rebsd(or(rebsd.prop.confidenceindex(:)<min_CI,rebsd.prop.confidenceindex(:)==-1)); %CI threshhold

    %% Plot and Save IPF Map
    fprintf('Plotting orientation data for %s...\n',fnameOut{num});  %Display name of current EBSD scan
    %Plot IPF Key
%     ipfKey = ipfHSVKey(rebsd);
%     figure; plot(ipfKey)
%     mtexTitle(' ');
%     print([pwd,'\Figures\','IPF Key'],'-dtiffn','-r350') ;
%     
    % Plot IPF Map
%     color = ipfKey.orientation2color(rebsd.orientations);
%     figure;plot(grains,'FaceColor','black','micronBar','off');
%     hold on
%     plot(rebsd,color,'micronBar','off')
%     plot(grains.boundary,'lineWidth',0.03,'linecolor','w','micronBar','off')
%     hold off
%     mtexTitle('IPF Map');
%     legend off
%     print([pwd,'\Figures\',fnameOut{num},' - IPF Map'],'-dtiffn','-r350');
%         
    %% Calc min+max grain axes (principle components - ellipse fit)
    fprintf('Calculating grain statistics for %s...\n',fnameOut{num});  %Display name of current EBSD scan
    %Subtract boundary grains
    grains_crop=grains;
    outerBoundary_id = any(grains_crop.boundary.grainId==0,2);  %find boundary grains
    grain_id = grains_crop.boundary(outerBoundary_id).grainId;  %compute grain_id for boundary grains
    grain_id(grain_id==0) = [];  %remove all zeroes
    grains_crop(grain_id) = [];  %remove boundary grains
    
    %Plot with micron bar
%     figure; plot(grains_crop)
%     mtexTitle('Grain Map w/o Boundary Grains');
%     legend off
%     print([pwd,'\Figures\',fnameOut{num},' - Micron Bar'],'-dtiffn','-r200') ;
%     
    %Calculate area weighted (AW) Grain diameters
    [omega,a,b]= principalComponents(grains_crop);
    a2 = a*2; %major diameter
    b2 = b*2; %minor diameter
    eqDiam = 2*equivalentRadius(grains_crop); %equivalent diameter
    
    a_avg = mean(a2);
    b_avg = mean(b2);
    eqDiam_avg = mean(eqDiam);
    
    AW = grains_crop.area./sum(grains_crop.area);
    repeat = round(AW/min(AW));
    dx = repelem(a2,repeat);
    dy = repelem(b2,repeat);
    dz = repelem(eqDiam,repeat);
    
    Dmaj = mean(dx);
    Dmin = mean(dy);
    Deq = mean(dz);
    
    
    %% Display
    fprintf('\nDavg_F_maj = %f um\n',a_avg);
    fprintf('Davg_F_min = %f um\n',b_avg);
    fprintf('Davg_AW_maj = %f um\n',Dmaj);
    fprintf('Davg_AW_min = %f um\n',Dmin);
    fprintf('Deq_AW = %f um\n\n',Deq);
    
    %% Save grain size data
    fprintf(fileID,[fnameOut{num},'\n']);
    fprintf(fileID,'==================================================\n');
    fprintf(fileID,'Mean grain MAJOR axis (FREQ. weight)  = %2.3f um\n',a_avg);
    fprintf(fileID,'Mean grain MINOR axis (FREQ. weight)  = %2.3f um\n',b_avg);
    fprintf(fileID,'Mean grain MAJOR axis (AREA weight)   = %2.3f um\n',Dmaj);
    fprintf(fileID,'Mean grain MINOR axis (AREA weight)   = %2.3f um\n',Dmin);
    fprintf(fileID,'Mean equiv. grain diam. (AREA weight) = %2.3f um\n',Deq);
    fprintf(fileID,'Scan size  = %16f, %16f um\n',max(rebsd.x),max(rebsd.y));
    fprintf(fileID,'==================================================\n\n');
    
    %% Options
    close all
    
end

fclose(fileID);

