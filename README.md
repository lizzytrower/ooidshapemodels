# ooidshapemodels
This suite of Matlab code accompanies the manuscript "Ooid cortical stratigraphy reveals the common histories of individual co-occuring sedimentary grains" by Trower et al.

Code for identifying cortical boundary layers:

ooidlampicking.m - aids user in identifying cortical layer boundaries from a thin section image

Code for modeling cortical boundary shapes corresponding to different growth-abrasion histories:

abrasioncalculator.m - calculates abrasion rate as a function of coordinates (theta,rho), Rouse number, and Stokes threshold

abrasioncalculatorH.m - same as abrasioncalculator, but including water depth (H) as additional input

fsimabrasion2.m - calculates new coordinates of grain cross section (theta,rho) after applying some abrasion increment

ooidgrowth.m - calculates new coordinates of grain cross section (theta,rho) after applying some growth increment

susp_abrasion_calculations_abrcalc.m - calculations used by abrasioncalculator.m

Code for plotting and visualizing results:

ooidhistvector.m - transforms model output of best-fit growth and abrasion increments into coordinates for a stair-step growth-abrasion history

ooidlamoverlay.m - plots best fit cortical layer boundaries over original image

plotooidhist.m - plots a stair-step plot of best-fit growth and abrasion steps

Example of growth-abrasion history simulation:

BB4cm_ROI13loop.m - code that uses the following two files to simulate growth-abrasion histories, then identify and plot the best-fit

BB4cm_ROI13data.mat - dataset with coordinates of cortical layer boundaries

BB4cm_ROI13.tif - corresponding thin section image

These codes require two other tools to run:

polygeom.m by H. J. Sommer, available on the Matlab File Exchange: https://www.mathworks.com/matlabcentral/fileexchange/319-polygeom-m

DiscreteFrechetDist.m by Zachary Danzier, available on the Matlab File Exchange: https://www.mathworks.com/matlabcentral/fileexchange/31922-discrete-frechet-distance
