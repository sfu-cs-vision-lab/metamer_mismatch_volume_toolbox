% Metamer Mismatch Volume Toolbox version 1.0
% This toolbox is available under GNU General Public version 3 license, https://www.gnu.org/licenses/gpl-3.0.en.html

% Any research publications which use any components of this toolbox in
% either original or modified form  should cite [1] and [2]:
% [1] Michal Mackiewicz, Hans Jakob Rivertz, and Graham Finlayson, 
%  "Spherical sampling methods for the calculation of metamer mismatch
%  volumes," J. Opt. Soc. Am. A 36, 96-104 (2019)
% [2] Graham D. Finlayson and Peter Morovic, "Metamer sets," J. Opt. Soc.
% Am. A 22, 810-819 (2005) 

% Michal Mackiewicz, University of East Anglia, 2021

% Calculates metamer mismatch volumes for grey-scale reflectance for a
% change of illumants.
% Should be read together with [1].
clear
resol = 380:1:735;
load data/T_xyzJuddVos
nw = numel(resol);
R = interp1(380:5:780,T_xyzJuddVos',resol,'cubic');

%Illuminant 1 D65;
%Illuminant 2 A;
load data/D65_380_1_735
E1 = interp1(380:1:735,E,resol,'cubic')';

L1 = diag(E1)*R;
load data/IllA
E2 = interp1(IllA(:,1),IllA(:,2),resol,'cubic')';
    
L2 = diag(E2)*R;
%%
refl_gr = .5*ones(size(E1))'; %.5 grey reflectance
ros = refl_gr*diag(E1)*R;
ros2 = refl_gr*diag(E2)*R;
ns = 10^6;   %number of samples for spherical sampling
rngset = rng;

ort_flag = 1;%[0 or 1]; see [1], Fig. 2 - orthonormal or standard sensors

%returns half-space representation of 6-D Object Colour Solid
[IneqCon,bIneqCon]=objectColSol_sphericalSampling([L1,L2],ns,rngset,ort_flag);

% intersection with equality constraint
[IneqCon,bIneqCon] = normalise_rows(IneqCon,bIneqCon);
EqCon = [1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0];

bEq = ros'; % we will calculate MMV for 0.5 grey
NullEqCon = null(EqCon);
    
x0=EqCon\bEq;
NewCon=IneqCon*NullEqCon;
bNew=bIneqCon-IneqCon*x0;

%calculates vertices of half-space intersection
vertices = calculateIntersectionVertices([NewCon,-bNew],ros2');
    
figure
[hull,v]=convhulln(vertices,{'Qt','C0.001'});
trisurf(hull,vertices(:,1),vertices(:,2),vertices(:,3));