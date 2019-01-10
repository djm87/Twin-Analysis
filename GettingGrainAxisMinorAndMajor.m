%**************************************************************************
% Determine all MTEX standard Grain Size and Shape parameters
%**************************************************************************
% calculate all parameters
% area (with holes = non-indexed points)
g_area = area(grains);
% aspect ratio is the ratio between the two principal components
g_aspectratio = aspectratio(grains);
% centroid = barycenter of the grain-polygon, with respect to its holes
% note that the barycenter can be outside the grain for some complex shapes
g_centroid = centroid(grains);
% equivalent radius er = sqrt(area/pi)
g_equivalentradius = equivalentradius(grains);
% equivalent diameter ed = 2*sqrt(area/pi)
g_equivalentdiameter = 2 * equivalentradius(grains);
% equivalent perimeter ep = 2*pi*equivalent radius
g_equivalentperimeter =  equivalentperimeter(grains);
% calculates the perimeter (p) of a grain, with holes
g_perimeter = perimeter(grains);          
% shapefactor = perimeter/equivalent perimeter = p/ep
g_shapefactor = shapefactor(grains);
%

There are two functions for grain shape; aspect ratio and shape factor, I prefer aspect ratio 
of course can make tests with both.  To get the orientation of the grain long axis associated with aspect ratio

There used to be a function for plotting the ellipse of the grain on EBSD map, but this was withdrawn from MTEX several years ago, it might be nice to have back.
%%
%**************************************************************************
% Determine angle of long axis of grain from EBSD map X-axis
%**************************************************************************
% Principle components using SVD analysis
% Eigen-vectors (ev) Eigen-values (ew)
[ev,ew]= principalcomponents(grains,'Hull');
% Angle from X-axis using Eigen-vectors (ev)
Theta = atan2(ev(1,2,:),ev(1,1,:))/degree;
% Remove unwanted singleton dimensions
Theta=squeeze(Theta);
% keep angles in the range 0-180 degrees from X-axis
for i=1:length(Theta)
    if(Theta(i)<0)
       Theta(i) = 180 + Theta(i);
    end
end
%%
%**************************************************************************
% Plot classical Rose diagram frequency vs angle
%**************************************************************************
figure
Rad_Theta = Theta * pi/180;
% bin width = 5 degrees (360/5 = 72)
n_angle_bins_Rose = 72;
rose(Rad_Theta,n_angle_bins_Rose)
% plot title
title('Angle from X=0 degrees Radius = counts','FontSize’,14)
% save plot
savefigure('/MatLab_Programs/Rose_diagram_angle_vs_aspect_ratio.png’)
        
        