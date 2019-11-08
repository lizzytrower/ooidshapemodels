function [thetaf,rhof,equivDf,IRf,areaf,Firey,drho1,dtheta,drho2,curvature] = ...
    fsimabrasion2( thetai,rhoi,abrasion_inc )
%FSIMABRASION2 Calculates change in grain shape due to abrasion
%   This function calculates the change in grain shape from some initial
%   outline (thetai, rhoi) to some final outline (thetaf, rhof), by
%   applying an abrasion increment (abrasion_inc).

%   (thetai,rhoi) is the initial outline, where thetai is in radians and
%   rhoi is in microns

%   abrasion_inc is a grain-averaged radial abrasion increment in microns

%   (thetaf,rhof) is the primary output - the outline after the abrasion
%   step, where thetaf is in radians and rhof is in microns

%   Additional outputs are available for troubleshooting errors. Note that
%   errors with fsimabrasion2 typically result from using an abrasion
%   increment that is too large.

%   This function was written by Lizzy Trower (University of Colorado
%   Boulder) in MATLAB 2018b on a Windows computer, last updated in
%   November 2019.

%set theta = [0,2pi]
negs = find(thetai < 0);
thetai(negs) = thetai(negs) + 2*pi;
dtheta = mean(diff(thetai));

%centered finite differences of first derivative
rhoi_back = cat(2,rhoi(end),rhoi(1:end-1));
rhoi_back2 = cat(2,rhoi(end-1:end),rhoi(1:end-2));
rhoi_fwd = cat(2,rhoi(2:end),rhoi(1));
rhoi_fwd2 = cat(2,rhoi(3:end),rhoi(1:2));
drho1 = (-rhoi_fwd2 + 8*rhoi_fwd - 8*rhoi_back + rhoi_back2)./(12*dtheta);

%centered finite differences of second derivative
drho2 = (-rhoi_fwd2 + 16*rhoi_fwd - 30*rhoi + 16*rhoi_back - rhoi_back2)./...
    (12*dtheta.^2);

curvature = (rhoi.^2 + 2.*drho1.^2 - rhoi.*drho2)./(rhoi.^2 + drho1.^2).^(3/2);

%average radial abrasion rate
A_ave = abrasion_inc; %equivalent radial abrasion increment (microns)

perimeter = sum((rhoi.^2+drho1.^2).^(1/2).*dtheta); %arc length from 0 to 2pi

Firey = 1 + perimeter.*(curvature);

normfactor = 2*pi*A_ave./sum(subplus(Firey).*dtheta);

rhonew = rhoi - subplus(Firey).*normfactor;

thetaf = thetai;
rhof = rhonew;

%set theta = [0,2pi]
negs1 = find(thetaf < 0);
thetaf(negs1) = thetaf(negs1) + 2*pi;

areaf = sum(dtheta.*rhof.^2./2);
equivDf = 2*sqrt(areaf/pi());
IRf = 4*pi()*areaf/perimeter^2;




