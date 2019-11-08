function [outputvector] = ooidhistvector(filename)
%OOIDHISTVECTOR Calculates ooid histories from model output
%   This function uses outputs from the ooid lamina simulation code
%   (specifically the best-fit growth increments and the # of iterations
%   for the best-fit abrasion surface) to plot a stair-step graph
%   illustrating cumulative growth (y axis) and abrasion time (x axis).
%   This enables one to graphically compare the growth-abrasion histories
%   of different ooids. Unlike 'plotooidhist', this code simply outputs the
%   (x,y) values rather than plotting them and is used for calculating the
%   Frechet distance matrix.

%   This function was written by Lizzy Trower (University of Colorado
%   Boulder) in MATLAB 2018b on a Windows computer, last updated in
%   November 2019.

load(filename,'growthinc_bestfits','ind_bestfits')

abrasionincs = ind_bestfits(:,2);
growthincs = growthinc_bestfits;
cgrowthincs = cumsum(growthincs);

deltat = 0.0001; %timestep in hr
abrasiontimes = abrasionincs*deltat; %[hr]
cabrasiontimes = cumsum(abrasiontimes);

%set up vectors for the (x,y) coordinates
xcoords = zeros(2*length(growthincs)+1,1);
ycoords = zeros(2*length(growthincs)+1,1);

for nn = 1:length(growthincs)
    ycoords(2*nn:2*nn+1) = cgrowthincs(nn);
    if nn == length(growthincs)
        xcoords(end) = cabrasiontimes(end);
    else
        xcoords(2*nn+1:2*nn+2) = cabrasiontimes(nn);
    end
end

outputvector = cat(2,xcoords,ycoords);

end

