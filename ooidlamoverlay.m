function [outputfig] = ooidlamoverlay(rhos_i,rhos_sim,...
    theta,imagename,xadj,yadj,calibratio)
%OOIDLAMOVERLAY Plots traced and simulated laminae on original image
%   This function plots the original lamina traces and simulated bestfit
%   traces over the original image, including re-scaling and translating
%   these traces so they match the image. Requires keeping track of the
%   xadj and yadj outputs from the ooidlampicking function.

%   This function was written by Lizzy Trower (University of Colorado
%   Boulder) in MATLAB 2018b on a Windows computer, last updated in
%   November 2019.

%rescale to plot on image, using original "calibratio" that scaled image
%traces into microns
rhos_i_recalib = rhos_i./calibratio;
rhos_sim_recalib = rhos_sim./calibratio;

%transform initial traces back into cartesian
xi = zeros(length(rhos_i(:,1)),length(rhos_i(1,:)));
yi = xi;
for ni = 1:length(rhos_i(:,1))
    [xtemp,ytemp] = pol2cart(theta,rhos_i_recalib(ni,:));
    xi(ni,:) = xtemp;
    yi(ni,:) = ytemp;
end

%translate coordinates to match image
xi = xi + xadj;
yi = yi + yadj;

%transform simulated traces into cartesian
xsim = zeros(length(rhos_sim(:,1)),length(rhos_sim(1,:)));
ysim = xsim;
for nsim = 1:length(rhos_sim(:,1))
    [xtemp2,ytemp2] = pol2cart(theta,rhos_sim_recalib(nsim,:));
    xsim(nsim,:) = xtemp2;
    ysim(nsim,:) = ytemp2;
end

%translate coordinates to match image;
xsim = xsim + xadj;
ysim = ysim + yadj;

%open image
im = imread(imagename);
outputfig = figure;
imshow(im)
hold on

%plot lamina traces in black
% for ni2 = 1:length(rhos_i(:,1))
%     plot(xi(ni2,:),yi(ni2,:),'k')
% end

%plot simulated lamina traces
for nsim2 = 1:length(rhos_sim(:,1))
    plot(xsim(nsim2,:),ysim(nsim2,:),'LineWidth',1)
end


end

