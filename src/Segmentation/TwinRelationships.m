% Based on A look-up table based approach to characterize
%crystal twinning for synchrotron X-ray Laue
%microdiffraction scans
%% Initialize
% clear all
coa=1.6250;

%% {10-11}<10-12> %contraction twin
%About {11-20}
acosd((-4*coa^2+3)/(4*coa^2+3))

%% {10-12}<10-11> %tensile twin
%About {11-20}
acosd((-coa^2+3)/(coa^2+3))

%% {10-13}<30-32>
%About {11-20}
acosd((-4*coa^2+27)/(4*coa^2+27))
acosd((-2*coa^2-1)/(2*coa^2+2))
%% {11-21}<11-26> Tensile Twin
%About {11-20}
acosd((-4*coa^2+1)/(4*coa^2+1))

%About {11-23}
acosd((2*coa^2-1)/(4*coa^2+1))

%About {10-10}
acosd((-4*coa^2+1)/(4*coa^2+1))
%% {11-22}<11-23> 
%About {10-10}
acosd((-coa^2+1)/(coa^2+1))

%About {22-43}
acosd((coa^2-2)/(2*coa^2+2))

%% {11-23}<11-22> present
%About {10-10}
acosd((-4*coa^2+9)/(4*coa^2+9))

%% {11-24}<22-43>
%About {10-10}
acosd((-coa^2+4)/(coa^2+4))