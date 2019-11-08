function [xcoords,ycoords] = plotooidhist(growthincs,abrasionincs,plotname)
%PLOTOOIDHIST Function to plot ooid history using simulation outputs
%   This function uses outputs from the ooid lamina simulation code
%   (specifically the best-fit growth increments and the # of iterations
%   for the best-fit abrasion surface) to plot a stair-step graph
%   illustrating cumulative growth (y axis) and abrasion time (x axis).
%   This enables one to graphically compare the growth-abrasion histories
%   of different ooids.

%   growthincs are the best-fit growth increments, in microns

%   abrasionincs are the best fit number of abrasion increments

%   plotname is a string with the name of the sample to add to the figure
%   legend

%   This function was written by Lizzy Trower (University of Colorado
%   Boulder) in MATLAB 2018b on a Windows computer, last updated in
%   November 2019.

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

gcf
hold on
plot(xcoords,ycoords)

h = findobj(gcf, 'Type', 'legend'); 
newstr = [h.String{1:end-1} {plotname}]; 
allDatah = flipud(get(gca,'children')); 
h.PlotChildren = allDatah; 
h.String = newstr;

end

