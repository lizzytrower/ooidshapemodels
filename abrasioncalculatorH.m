function [abrasionrate] = abrasioncalculatorH(thetain,rhoin,Rouse,Stokes,H)
%ABRASIONCALCULATORH Calculates radial abrasion rate for an ooid
%   This function calculates radial abrasion rate for an ooid as a function
%   of the ooid's size (defined in polar coordinates), Rouse number, Stokes
%   threshold, and water depth (H)

%   (thetain,rhoin) are the coordinates that describe the 2D cross-section
%   of the grain, where thetain is in radians and rhoin is in microns

%   H is water depth in meters

%   the output, 'abrasionrate', is in units of microns/hour

%   This function was written by Lizzy Trower (University of Colorado
%   Boulder) in MATLAB 2018b on a Windows computer, last updated in
%   November 2019.

[x0,y0] = pol2cart(thetain,rhoin);

[geo_i] = polygeom(x0,y0);
equivD = 2*sqrt(geo_i(1)/pi())*10^-6; %[m]
rho_s = 2800; %[kg/m^3] density of aragonite
rho_f = 1080; %[kg/m^3] density of Great Salt Lake water
R=(rho_s - rho_f)/rho_f; %[dimensionless] submerged specific density
g = 9.81; %[m/s^2]
young = 144*10^9; %[Pa] young's modulus
strength = 1*10^6; %[Pa] tensile strength
kv = 400*10^4; %[dimensionless]
tauc = 0.03; %Critical Shields number.  0.03 is good for sand.
Stc = Stokes; %critical Stokes number
nu = 1.3*10^-6; %[m^2/s] kinematic viscosity of water
CSF = 0.8;  %1 is for spheres, 0.8 is for natural
PS = 3.5;  %6 is for spheres, 3.5 is for natural
Dstar = (R.*g.*equivD.^3)./(nu.^2);
X = log10(Dstar);
R1 = -3.76715+1.92944.*X - 0.09815.*(X.^2) - 0.00575.*(X.^3) + 0.00056.*(X.^4);
R2 = log10(1-((1-CSF)./0.85))-(((1-CSF).^2.3).*tanh(X-4.6)) + 0.3.*(0.5-CSF).*((1-CSF).^2).*(X-4.6);
R3 = (0.65-((CSF./2.83).*tanh(X-4.6))).^(1+((3.5-PS)./2.5));
Wstar = R3.*10.^(R2+R1);
ws = (R.*g.*nu.*Wstar).^(1./3); %[m/s]
cdrag = (4/3).*(R.*g.*equivD)./(ws.^2);
betta = 1;
ustar = ws./(Rouse.*.41.*betta); %[m/s]
D = equivD; %[m]
tau = ustar^2/(R*g*D); %[dimensionless]
tstage = tau/tauc; %[dimensionless]
A1 = 0.36; %[dimensionless]
susp_abrasion_calculations_abrcalc
abrasionrate = Rabrasion; %[um/hr]
end

