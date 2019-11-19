%   This code is an example workflow to simulate possible ooid cortical
%   layer boundaries resulting from different growth-abrasion histories,
%   identify the best-fit history, and overlay it on the original thin
%   section image. The code requires an input file (here,
%   BB4cm_ROI13data.mat) with vectors of traced cortical layer boundaries
%   (in polar coordinates) and the (x,y) coordinates of the grain's center
%   in the original input image (for the overlay plot). Because abrasion
%   requires very small time steps (and is therefore time-consuming to
%   run), the following code is designed to use Matlab's parallel computing
%   toolbox.

%   This code was written by Lizzy Trower (University of Colorado Boulder)
%   in Matlab 2018b on a Windows computer, last updated in November 2019.


clear
tic
load("BB4cm_ROI13data.mat")

deltat = 0.0001; %timestep in hr
Rouse = 2.5; %Rouse number
Stc = 1; %critical Stokes number

repeat = 400000; %number of abrasion steps

%create array to store outputs in
rho_bestfits = zeros(length(BB4cm_ROI13.lams(:,1))-1,(length(BB4cm_ROI13.lams(1,:))-1),'single');
growth_bestfits = zeros(length(BB4cm_ROI13.lams(:,1))-1,(length(BB4cm_ROI13.lams(1,:))-1),'single');
growthinc_bestfits = zeros(length(BB4cm_ROI13.lams(:,1))-1,1,'single');
ind_bestfits = zeros(length(BB4cm_ROI13.lams(:,1))-1,2,'single');

for nnn = 1:length(BB4cm_ROI13.lams(:,1))-1
    
    theta0 = single(BB4cm_ROI13.theta);
    rho0 = zeros(1,length(theta0),'single');
    
    if nnn == 1
        rho0 = BB4cm_ROI13.lams(nnn,:);
    else
        rho0(1:end-1) = rho_bestfits(nnn-1,:);
        rho0(end) = rho0(end-1);
    end
    rho_matchlam = single(BB4cm_ROI13.lams(nnn+1,:));
    
    mininc = round(min(rho_matchlam - rho0));
    growthinc = mininc:5:mininc+20; %[um] growth increments

    [thetas_growth] = zeros(length(growthinc),length(theta0),'single');
    [rhos_growth] = zeros(length(growthinc),length(rho0),'single');
    rhos_abrade = zeros(length(growthinc),length(rho0)-1,repeat,'single');
    misfit_sqsum = zeros(length(growthinc),repeat,'single');
    misfit_goal = rho_matchlam(1:length(theta0) - 1);
    
    parfor ii = 1:length(growthinc)
    
    [thetaii,rhoii] = ooidgrowth(theta0,rho0,growthinc(ii),1);
    replaceval = find(isfinite(rhoii));
    rhoii(isnan(rhoii)) = rhoii(replaceval(1));
    thetas_growth(ii,:) = thetaii';
    rhos_growth(ii,:) = rhoii';
    
    end
    
    parfor aa = 1:length(growthinc)

    for jj = 1:repeat
        
        k = jj/1000;
        
        if jj == 1   
            saveind = 0;
            theta_working = thetas_growth(aa,1:length(theta0)-1); %#ok<PFBNS>
            rho_working = rhos_growth(aa,1:length(theta0)-1); %#ok<PFBNS>       
        end

        if jj == 1
            abrasionrate = abrasioncalculator(theta_working,rho_working,...
                Rouse,Stc);
        elseif floor(k) == k
            abrasionrate = abrasioncalculator(theta_working,rho_working,...
                Rouse,Stc);
        end
        
        abrasion_inc = abrasionrate.*deltat;
        
        [theta_working,rho_working] = ...
            fsimabrasion2(theta_working,rho_working,abrasion_inc);

        rhos_abrade(aa,:,jj) = rho_working;
        misfits = rho_working - misfit_goal; 
        misfits_sq = misfits.^2;
        misfit_sqsum(aa,jj) = sum(misfits_sq);
        
    end
    
    end
    
    [bfval,bfind] = min(misfit_sqsum(:));
    [bf_row, bf_col] = ind2sub(size(misfit_sqsum),bfind);
    rho_bestfits(nnn,:) = rhos_abrade(bf_row,:,bf_col);
    growth_bestfits(nnn) = rhos_growth(bf_row,1:end-1);
    ind_bestfits(nnn,1) = bf_row;
    ind_bestfits(nnn,2) = bf_col;
    
end
toc

save('BB4cm_ROI13_simdata_senstest.mat','ind_bestfits','rho_bestfits','growth_bestfits',...
    'growthinc_bestfits');
ooidlamoverlay(BB4cm_ROI13.lams(:,1:end-1),rho_bestfits,...
    BB4cm_ROI13.theta(1:end-1),'BB4cm_ROI13.tif',BB4cm_ROI13.xadj,BB4cm_ROI13.yadj,...
    .225)
saveas(gcf,'BB4cm_ROI13overlay.fig')
close
