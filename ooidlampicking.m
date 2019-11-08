function [ooidlampoints,theta_even,xadj,yadj] = ooidlampicking(imagename,numoflams,...
    calibratio)
%OOIDLAMPICKING Function to pick out lamina traces in ooids
%   Function that opens the selected image and allows user to pick out the
%   specified number of lamina traces in the image, which are then
%   converted to evenly-spaced polar coordinates. Number of lamina to trace
%   must be specified ("numoflams"); the scale of the image is adjusted
%   using "calibratio", which can be calculated using the following code:

% imfile = imagename;
% im = imread(imfile);
% fig = figure;
% imshow(im);
% calibmeas = imdistline;
% caliblength = getDistance(calibmeas);
% calibratio = calibline/caliblength;

%   This function was written by Lizzy Trower (University of Colorado
%   Boulder) in MATLAB 2018b on a Windows computer, last updated in
%   November 2019.

%open image
imfile = imagename;
im = imread(imfile);
fig = figure;
imshow(im);

%use @drawassisted to trace surface, then pick the x,y coordinates
lamnuc = drawassisted;
xptsnuc = lamnuc.Position(:,1);
yptsnuc = lamnuc.Position(:,2);

%finds centroids of initial shape
ooid_initial_geom = polygeom(xptsnuc,yptsnuc);
xadj = ooid_initial_geom(2);
yadj = ooid_initial_geom(3);

%translates shapes so centroid of initial is at (0,0)
xptsnuc_adj = xptsnuc - xadj;
yptsnuc_adj = yptsnuc - yadj;

%adjusts theta range to [0, 2pi]
[thetanuc,rhonuc] = cart2pol(xptsnuc_adj,yptsnuc_adj);
negsnuc = find(thetanuc < 0);
thetanuc(negsnuc) = thetanuc(negsnuc) + 2*pi;

%sorts (theta,rho) pairs so theta is monotonically increasing
[minnuc,startnuc] = min(thetanuc);
thetanuc_adj = cat(1,thetanuc(startnuc:end),thetanuc(1:startnuc-1));
rhonuc_adj = cat(1,rhonuc(startnuc:end),rhonuc(1:startnuc-1));

%removes duplicate points
[vals_adj] = unique([rhonuc_adj thetanuc_adj],'stable','rows');

%create evenly spaced trace w/ 2000 points
theta_even = linspace(0,2*pi,2000);
rhonuc_even = interp1(vals_adj(:,2),vals_adj(:,1),theta_even);

%find and replace NaNs
replaceval_nuc = find(isfinite(rhonuc_even));
rhonuc_even(isnan(rhonuc_even)) = rhonuc_even(replaceval_nuc(1));

%create matrix to add additional lamina traces to
ooidlampoints = zeros(numoflams,length(theta_even));
ooidlampoints(1,:) = rhonuc_even;

%forloop to trace and transform each lamina for the selected number of
%iterations ("numoflams")
for n = 2:numoflams
    %get the trace of the lamina
    lamtrace = drawassisted;
    xpts = lamtrace.Position(:,1);
    ypts = lamtrace.Position(:,2);
    
    %adjust so centered @ origin, w/r/t nucleus centroid
    xpts_adj = xpts - ooid_initial_geom(2);
    ypts_adj = ypts - ooid_initial_geom(3);
    
    %transform to polar coordinates and sort so theta is monotonically
    %increasing
    [theta,rho] = cart2pol(xpts_adj,ypts_adj);
    negs = find(theta < 0);
    theta(negs) = theta(negs) + 2*pi;
    [minval,startpt] = min(theta);
    theta_adj = cat(1,theta(startpt:end),theta(1:startpt-1));
    rho_adj = cat(1,rho(startpt:end),rho(1:startpt-1));
    
    %remove duplicate points
    [dups_adj] = unique([rho_adj theta_adj],'stable','rows');
    
    %create evenly spaced trace w/ 2000 points (same theta values as for
    %nucleus)
    rho_even = interp1(dups_adj(:,2),dups_adj(:,1),theta_even); 
    
    %find and replace NaNs
    replaceval = find(isfinite(rho_even));
    rho_even(isnan(rho_even)) = rho_even(replaceval(1));
    ooidlampoints(n,:) = rho_even;
end

%resize using the "calibratio" input
ooidlampoints = ooidlampoints*calibratio;

end

