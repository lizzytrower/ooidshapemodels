function [ thetaf,rhof ] = ooidgrowth( thetai,rhoi,increment,nrepeats )
%OOIDGROWTH Function that applies 2D surface-normal growth
%   This function takes a set of polar coordinates that define a 2D shape
%   (thetai, rhoi) and adds a certain number (nrepeats) of surface-normal
%   layers.

%   This function was written by Lizzy Trower (University of Colorado
%   Boulder) in MATLAB 2018b on a Windows computer, last updated in
%   November 2019.

%the math is easier if we convert back into Cartesian (x,y) coordinates
[xi,yi] = pol2cart(thetai,rhoi);

%this code works by taking each point along the initial shape and drawing a
%circle with the radius of the growth increment, so here we define that
%radius and set up a set of angles over which to draw that circle
th_int = 0:0.1:2*pi;
radius = increment;

%optional code to plot a figure of the initial shape:
% figure
% polarplot(thetai,rhoi)
% hold on

for m = 1:nrepeats
    
    if m == 1
        xi = xi;
        yi = yi;
    else
        xi = xf_ch;
        yi = yf_ch;
    end

    xx = ones(length(xi),length(th_int));
    yy = ones(length(xi),length(th_int));

    %here is where we draw all the circles
    for n = 1:length(xi)
        xx(n,:) = radius.*cos(th_int) + xi(n);
        yy(n,:) = radius.*sin(th_int) + yi(n);
    end

    %combine all the points from the circles into single vectors
    xf = xx(:);
    yf = yy(:);

    %this function finds the outer boundary of all the circles that we just
    %drew, which becomes our new shape
    kf = boundary(xf,yf,.4);

    %these will get used as inputs for the next time step if nrepeats>1
    xf_ch = xf(kf);
    yf_ch = yf(kf);

    [thetaf,rhof] = cart2pol(xf_ch,yf_ch);
    %optional code to plot this new surface:
%     polarplot(thetaf,rhof)
    
    %everything works better if we only work at angles between 0 and 2pi,
    %so these lines fix that
    negs = find(thetaf < 0);
    thetaf(negs) = thetaf(negs) + 2*pi;

    %this section adds a fix in case there are any weird loops in our new
    %shape
    spacing = diff(thetaf);
    th_nonzero = find(spacing);
    thetaf = cat(1,thetaf(th_nonzero),thetaf(end));
    rhof = cat(1,rhof(th_nonzero),rhof(end));
    
    theta_even = linspace(0,2*pi,length(thetai));
    rho_even = interp1(thetaf(1:end-1),rhof(1:end-1),theta_even);
    
    thetaf = theta_even;
    rhof = rho_even;
 
end

end