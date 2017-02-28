function npptotalmax = hpm_vegNPP20(wtd_opt,wtd_range,pd_opt,pd_range,npp_rel,k0,num_veg)
% function tot_npp = hpm_vegNPP12(wtd_opt,wtd_range,pd_opt,pd_range,npp_rel,k0)

% v.20: making version for version control - June 2015

% v12 for arbitary number of PFTs

% v9

% v8 no change from v6

% function generates plots of veg NPP and calculates total NPP

% each vegetation type will have an optimal position, and a surface of diminishing NPP

% response modeled as 2-D Gaussian functions, with (in some cases)
% different variances on either side of the optimum to skew the function
% to have a long (productivity) tail in one direction (e.g., increasing WTD)
% and a short (productivity) tail in the other direction (e.g., decreasing WTD)

% PARAMETERS
%-----------

%  WTD_opt = optimum WTD
%  a1 = NPP sensitivity to WTD increasing (deeper WT) (Gaussian curve variance)
%  a2 = NPP sensitivity to WTD decreasing (shallower WT) (Gaussian curve variance)

%  PD_opt = optimum PD or ALD
%  b1 = NPP sensitivity to PD increasing (deeper peat profile) (Gaussian curve variance)
%  b2 = NPP sensitivity to PD decreasing (Gaussian curve variance)

%  c = relative maximum NPP

[X,Y] = meshgrid(-.5:0.1:1, 0:0.1:8);
minarray = 0.00001*ones(size(X));
NPPmax = zeros(13,1);

XX = [-0.5:0.1:1];
YY = [0:0.1:8];

% loop through PFTs

for nveg = 1:1:num_veg

    atest1 = XX < wtd_opt(nveg);  % one if WTD <  optimum, zero otherwise
    atest2 = XX >= wtd_opt(nveg); % one if WTD >= optimum, zero otherwise
    btest1 = YY < pd_opt(nveg);   % one if PD  <  optimum, zero otherwise
    btest2 = YY >= pd_opt(nveg);  % one if PD  >= optimum, zero otherwise
    btest1 = btest1';
    btest2 = btest2';
    
    NPPa1 = exp(-(((XX - wtd_opt(nveg)) / wtd_range(1,nveg)).^2));
    NPPa2 = exp(-(((XX - wtd_opt(nveg)) / wtd_range(2,nveg)).^2));
    NPPb1 = exp(-((YY - pd_opt(nveg)) / pd_range(1,nveg)).^2);
    NPPb2 = exp(-((YY - pd_opt(nveg)) / pd_range(2,nveg)).^2);
    NPPb1 = NPPb1';
    NPPb2 = NPPb2';
    
    NPP1 = npp_rel(nveg) * ((NPPb1 .* btest1 + NPPb2 .* btest2) * (NPPa1 .* atest1 + NPPa2 .* atest2));
    NPPmax(nveg) = max(max(NPP1));
    NPP(nveg,:,:) = max(NPP1, minarray);
    
    if (nveg == 1)
        NPPtotal = NPP(1,:,:);
    else
        NPPtotal = NPPtotal + NPP(nveg,:,:);
    end
end

% TOTAL OF PFTs 1 to 12 (NO TREES); TOTAL OF MOSSES; TOTAL OF VASCULAR (NO
% TREES) 


% NPPtotal = NPP(1,:,:) + NPP(2,:,:) + NPP(3,:,:) + NPP(4,:,:) + NPP(5,:,:) + NPP(6,:,:) + ...
%            NPP(7,:,:) + NPP(8,:,:) + NPP(9,:,:) + NPP(10,:,:) + NPP(11,:,:) + NPP(12,:,:) + NPP(13,:,:);
% NPPmosstotal = NPP(5,:,:) + NPP(6,:,:) + NPP(7,:,:) + NPP(8,:,:) + NPP(9,:,:);
% NPPvasctotal = NPP(1,:,:) + NPP(2,:,:) + NPP(3,:,:) + NPP(4,:,:) +  ...
%            NPP(10,:,:) + NPP(11,:,:) + NPP(12,:,:); + NPP(13,:,:);
% NPPminerototal = NPP(1,:,:) + NPP(2,:,:) + NPP(4,:,:) + NPP(5,:,:) +  ...
%            NPP(3,:,:);
% NPPombrototal = NPP(6,:,:) + NPP(7,:,:) + NPP(8,:,:) + NPP(9,:,:) +  ...
%            NPP(10,:,:) + NPP(11,:,:) + NPP(12,:,:);  
% NPPtreestotal = NPP(13,:,:);      
% NPPminerofrac = NPPminerototal ./(NPPminerototal + NPPombrototal + eps);
% 
% NPPpeattotal = NPP(1,:,:)/k0(1) + NPP(2,:,:)/k0(2) + NPP(3,:,:)/k0(3) + ...
%                NPP(4,:,:)/k0(4) + NPP(5,:,:)/k0(5) + NPP(6,:,:)/k0(6) + ...
%                NPP(7,:,:)/k0(7) + NPP(8,:,:)/k0(8) + NPP(9,:,:)/k0(9) + ...
%                NPP(10,:,:)/k0(10) + NPP(11,:,:)/k0(11) + ...
%                NPP(12,:,:)/k0(12)  + NPP(13,:,:)/k0(13);
% 
% NPPgrassfrac = NPP(1,:,:) ./ (NPPtotal + eps);
% NPPminforbfrac = NPP(2,:,:) ./ (NPPtotal + eps);
% NPPminsedgefrac = NPP(3,:,:) ./ (NPPtotal + eps);
% NPPminshrubfrac = NPP(4,:,:) ./ (NPPtotal + eps);
% NPPbrownmossfrac = NPP(5,:,:) ./ (NPPtotal + eps);
% NPPholsphagfrac = NPP(6,:,:) ./ (NPPtotal + eps);
% NPPlawnsphagfrac = NPP(7,:,:) ./ (NPPtotal + eps);
% NPPhumsphagfrac = NPP(8,:,:) ./ (NPPtotal + eps);
% NPPfeatherfrac = NPP(9,:,:) ./ (NPPtotal + eps);
% NPPombforbfrac = NPP(10,:,:) ./ (NPPtotal + eps);
% NPPombsedgefrac = NPP(11,:,:) ./ (NPPtotal + eps);
% NPPombshrubfrac = NPP(12,:,:) ./ (NPPtotal + eps);
% NPPtreesfrac = NPP(13,:,:) ./ (NPPtotal + eps);
% 
NPPtotal = squeeze(NPPtotal);
% NPPpeattotal = squeeze(NPPpeattotal);
% NPPmosstotal = squeeze(NPPmosstotal);
% NPPvasctotal = squeeze(NPPvasctotal);
% NPPminerototal = squeeze(NPPminerototal);
% NPPombrototal = squeeze(NPPombrototal);
% NPPtreestotal = squeeze (NPPtreestotal);
% NPPminerofrac = squeeze(NPPminerofrac);
% NPPtotal = max(NPPtotal, minarray);
% NPPmosstotal = max(NPPmosstotal, minarray);
% NPPvasctotal = max(NPPvasctotal, minarray);
% NPPminerototal = max(NPPminerototal, minarray);
% NPPombrototal = max(NPPombrototal, minarray);
% NPPtreestotal = max (NPPtreestotal, minarray);

% NPPgrassfrac = squeeze(NPPgrassfrac);
% NPPminforbfrac = squeeze(NPPminforbfrac);
% NPPminsedgefrac = squeeze(NPPminsedgefrac);
% NPPminshrubfrac = squeeze(NPPminshrubfrac);
% NPPbrownmossfrac = squeeze(NPPbrownmossfrac);
% NPPholsphagfrac = squeeze(NPPholsphagfrac);
% NPPlawnsphagfrac = squeeze(NPPlawnsphagfrac);
% NPPhumsphagfrac = squeeze(NPPhumsphagfrac);
% NPPfeatherfrac = squeeze(NPPfeatherfrac);
% NPPombforbfrac = squeeze(NPPombforbfrac);
% NPPombsedgefrac = squeeze(NPPombsedgefrac);
% NPPombshrubfrac = squeeze(NPPombshrubfrac);
% NPPtreesfrac = squeeze(NPPtreesfrac); 
% 
npptotalmax = max(max(NPPtotal))
% npppeattotalmax = max(max(NPPpeattotal));
% nppmossmax = max(max(NPPmosstotal))
% nppvascmax = max(max(NPPvasctotal))
% nppmineromax = max(max(NPPminerototal))
% nppombromax = max(max(NPPombrototal))
% npptreesmax = max(max(NPPtreestotal))

% plot results

% plot or not

plot_flag = 1;

if (plot_flag > 0)

  % ---------  
  % figure(21)
  % ---------  

    figure(21)
    
    for nveg = 1:1:num_veg
        subplot(4,4,nveg)
        contourf(X, Y, squeeze(NPP(nveg,:,:)),100,'LineStyle','none')
        xlim([-0.25 1])

%         str = sprintf('PFT #', nveg);
%         title(str, 'FontSize',14)   
        title(['\fontsize{14}PFT #',num2str(nveg)]);

        ylim([0 8])
        zlim([0 1])
        colorbar
        caxis([0 2])
        xlabel('\fontsize{14}water table depth [m]')
        ylabel('\fontsize{14}peat depth [m]')
        zlabel('\fontsize{14}relative NPP')
        h1a=gca;
        set(h1a,'FontSize',14)
    end
    
    subplot(4,4,nveg+1)
    contourf(X, Y, NPPtotal,100,'LineStyle','none')
    xlim([-0.25 1])
    title('\fontsize{14}all PFT total');
    ylim([0 8])
    zlim([0 1])
    colorbar
    caxis([0 2])
    xlabel('\fontsize{14}water table depth [m]')
    ylabel('\fontsize{14}peat depth [m]')
    zlabel('\fontsize{14}relative NPP')
    h1a=gca; 
    set(h1a,'FontSize',14)

end

return

% **************************************************
% *** THE REST OF THIS IS OLD STUFF AND NOT USED ***    
% **************************************************

%     subplot(4,4,2)
%     contourf(X, Y, squeeze(NPP(2,:,:)),100,'LineStyle','none')
%     title('\fontsize{14}minerotrophic forb');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 2])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1b=gca; 
%     set(h1b,'FontSize',14)
% 
%     subplot(4,4,3)
%     contourf(X, Y, squeeze(NPP(3,:,:)),100,'LineStyle','none')
%     title('\fontsize{14}minerotrophic sedge');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 2])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
% 
%     subplot(4,4,4)
%     contourf(X, Y, squeeze(NPP(4,:,:)),100,'LineStyle','none')
%     title('\fontsize{14}minerotrophic shrub');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 2])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1a=gca; 
%     set(h1a,'FontSize',14)
%     
%     subplot(4,4,5)
%     contourf(X, Y, squeeze(NPP(5,:,:)),100,'LineStyle','none')
%     title('\fontsize{14}brown moss');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 2])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1b=gca; 
%     set(h1b,'FontSize',14)
%     
%     subplot(4,4,6)
%     contourf(X, Y, squeeze(NPP(6,:,:)),100,'LineStyle','none')
%     title('\fontsize{14}hollow sphagnum');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 1])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
% 
%     subplot(4,4,7)
%     contourf(X, Y, squeeze(NPP(7,:,:)),100,'LineStyle','none')
%     title('\fontsize{14}lawn sphagnum');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 1])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
% 
%     subplot(4,4,8)
%     contourf(X, Y, squeeze(NPP(8,:,:)),100,'LineStyle','none')
%     title('\fontsize{14}hummock sphagnum');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 1])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
%     
%     subplot(4,4,9)
%     contourf(X, Y, squeeze(NPP(9,:,:)),100,'LineStyle','none')
%     title('\fontsize{14}feather moss');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 1])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
%     
%     subplot(4,4,10)
%     contourf(X, Y, squeeze(NPP(10,:,:)),100,'LineStyle','none')
%     title('\fontsize{14}ombrotrophic forb');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 1])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
%     
%     subplot(4,4,11)
%     contourf(X, Y, squeeze(NPP(11,:,:)),100,'LineStyle','none')
%     title('\fontsize{14}ombrotrophic sedge');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 1])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
%     
%     subplot(4,4,12)
%     contourf(X, Y, squeeze(NPP(12,:,:)),100,'LineStyle','none')
%     title('\fontsize{14}ombrotrophic shrub');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 1])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
%     
%     subplot(4,4,13)
%     contourf(X, Y, squeeze(NPP(13,:,:)),100,'LineStyle','none')
%     title('\fontsize{14}trees');
%     xlim([-0.5 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 1])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)

% end 

% plot_flag = -1;


%   % ---------  
%   % figure(22)
%   % ---------  
% 
%   if (plot_flag > 0)
%       
%     figure(22)
%     subplot(4,1,1)
%     contourf(X, Y, NPPtotal,100,'LineStyle','none')
%     title('\fontsize{14}all PFTs');
%     xlim([-0.25 1])
%     ylim([0 8])
%     %zlim([0 5])
%     colorbar
%     caxis([0 6])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     %zlabel('\fontsize{14}relative NPP')
%     h2a=gca; 
%     set(h2a,'FontSize',14)
%     
%     subplot(4,1,2)
%     contourf(X, Y, NPPpeattotal,100,'LineStyle','none')
%     title('\fontsize{14}all PFTs - relative PEAT potential production');
%     xlim([-0.25 1])
%     ylim([0 8])
%     %zlim([0 5])
%     colorbar
%     caxis([0 npppeattotalmax])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     %zlabel('\fontsize{14}relative PEAT NPP')
%     h2a=gca; 
%     set(h2a,'FontSize',14)
%     
%     subplot(4,1,3)
%     contourf(X, Y, NPPminerofrac,100,'LineStyle','none')
%     title('\fontsize{14}minerotrophic fraction of total NPP');
%     xlim([-0.25 1])
%     ylim([0 8])
%     %zlim([0 5])
%     colorbar
%     caxis([0 1])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     %zlabel('\fontsize{14}minerotrophic fraction of NPP')
%     h2a=gca; 
%     set(h2a,'FontSize',14)
%     
%     subplot(4,1,4)
%     contourf(X, Y, NPPtreesfrac,100,'LineStyle','none')
%     title('\fontsize{14}trees fraction of total NPP');
%     xlim([-0.25 1])
%     ylim([0 8])
%     %zlim([0 5])
%     colorbar
%     caxis([0 1])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     %zlabel('\fontsize{14}minerotrophic fraction of NPP')
%     h2a=gca; 
%     set(h2a,'FontSize',14)
%     
%   end
%   
%   % ---------  
%   % figure(23)
%   % ---------  
% 
%   if (plot_flag > 0)
% 
%     figure(23)
%     subplot(3,2,1)
%     contourf(X, Y, NPPminerototal,100,'LineStyle','none')
%     title('\fontsize{14}minerotrophic PFTs');
%     xlim([-0.25 1])
%     ylim([0 8])
%     %zlim([0 5])
%     colorbar
%     caxis([0 6])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     %zlabel('\fontsize{14}relative NPP')
%     h2a=gca; 
%     set(h2a,'FontSize',14)
%     
%     subplot(3,2,2)
%     contourf(X, Y, NPPombrototal,100,'LineStyle','none')
%     title('\fontsize{14}ombrotrophic PFTs');
%     xlim([-0.25 1])
%     ylim([0 8])
%     %zlim([0 5])
%     colorbar
%     caxis([0 3])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     %zlabel('\fontsize{14}relative NPP')
%     h2a=gca; 
%     set(h2a,'FontSize',14)
%     
%     subplot(3,2,3)
%     contourf(X, Y, NPPmosstotal,100,'LineStyle','none')
%     title('\fontsize{14}all moss/lichen PFTs');
%     xlim([-0.25 1])
%     ylim([0 8])
%     %zlim([0 5])
%     colorbar
%     caxis([0 2])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     %zlabel('\fontsize{14}relative NPP')
%     h2a=gca; 
%     set(h2a,'FontSize',14)
%     
%     subplot(3,2,4)
%     contourf(X, Y, NPPvasctotal,100,'LineStyle','none')
%     title('\fontsize{14}all vascular PFTs');
%     xlim([-0.25 1])
%     ylim([0 8])
%     %zlim([0 5])
%     colorbar
%     caxis([0 6])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     %zlabel('\fontsize{14}relative NPP')
%     h2a=gca; 
%     set(h2a,'FontSize',14)
%     
%     subplot(3,2,5)
%     contourf(X, Y, NPPtreestotal,100,'LineStyle','none')
%     title('\fontsize{14}trees');
%     xlim([-0.25 1])
%     ylim([0 8])
%     %zlim([0 5])
%     colorbar
%     caxis([0 6])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     %zlabel('\fontsize{14}relative NPP')
%     h2a=gca; 
%     set(h2a,'FontSize',14)
%     
%   end
%   
% %  plot_flag = -1;
%   
%    % ---------  
%   % figure(24)
%   % ---------  
% 
%   if (plot_flag > 0)
% 
% 
%     figure(24)
%     subplot(4,1,1)
%     contourf(X, Y, squeeze(NPP(3,:,:)),100,'LineStyle','none')
%     title('\fontsize{14}minerotrophic sedge');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 2])
%   %  xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
% %    zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
%     
%     subplot(4,1,2)
%     contourf(X, Y, squeeze(NPP(8,:,:)),100,'LineStyle','none')
%     title('\fontsize{14}hummock sphagnum');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 1])
%  %   xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
% %    zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
%     
%     subplot(4,1,3)
%     contourf(X, Y, squeeze(NPP(12,:,:)),100,'LineStyle','none')
%     title('\fontsize{14}ombrotrophic shrub');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 1])
% %    xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
% %    zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
%     
%     subplot(4,1,4)
%     contourf(X, Y, NPPvasctotal,100,'LineStyle','none')
%     title('\fontsize{14}all vascular PFTs');
%     xlim([-0.25 1])
%     ylim([0 8])
%     %zlim([0 5])
%     colorbar
%     caxis([0 6])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     %zlabel('\fontsize{14}relative NPP')
%     h2a=gca; 
%     set(h2a,'FontSize',14)
%     
%   end
% 
% % plot or not
% 
% % plot_flag = 0;
% 
% if (plot_flag > 0)
% 
%   % ---------  
%   % figure(25)
%   % ---------  
% 
%     colorscale = 0.7;
%     
%     figure(25)
%     subplot(4,3,1)
%     contourf(X, Y, NPPgrassfrac,100,'LineStyle','none')
%     xlim([-0.25 1])
%     title('\fontsize{14}grass');
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 colorscale])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1a=gca; 
%     set(h1a,'FontSize',14)
% 
%     subplot(4,3,2)
%     contourf(X, Y, NPPminforbfrac,100,'LineStyle','none')
%     title('\fontsize{14}minerotrophic forb');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 colorscale])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1b=gca; 
%     set(h1b,'FontSize',14)
% 
%     subplot(4,3,3)
%     contourf(X, Y, NPPminsedgefrac,100,'LineStyle','none')
%     title('\fontsize{14}minerotrophic sedge');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 colorscale])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
% 
%     subplot(4,3,4)
%     contourf(X, Y, NPPminshrubfrac,100,'LineStyle','none')
%     title('\fontsize{14}minerotrophic shrub');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 colorscale])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1a=gca; 
%     set(h1a,'FontSize',14)
%     
%     subplot(4,3,5)
%     contourf(X, Y, NPPbrownmossfrac,100,'LineStyle','none')
%     title('\fontsize{14}brown moss');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 colorscale])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1b=gca; 
%     set(h1b,'FontSize',14)
%     
%     subplot(4,3,6)
%     contourf(X, Y, NPPholsphagfrac,100,'LineStyle','none')
%     title('\fontsize{14}hollow sphagnum');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 colorscale])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
% 
%     subplot(4,3,7)
%     contourf(X, Y, NPPlawnsphagfrac,100,'LineStyle','none')
%     title('\fontsize{14}lawn sphagnum');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 colorscale])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
% 
%     subplot(4,3,8)
%     contourf(X, Y, NPPhumsphagfrac,100,'LineStyle','none')
%     title('\fontsize{14}hummock sphagnum');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 colorscale])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
%     
%     subplot(4,3,9)
%     contourf(X, Y, NPPfeatherfrac,100,'LineStyle','none')
%     title('\fontsize{14}feather moss');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 colorscale])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
%     
%     subplot(4,3,10)
%     contourf(X, Y, NPPombforbfrac,100,'LineStyle','none')
%     title('\fontsize{14}ombrotrophic forb');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 colorscale])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
%     
%     subplot(4,3,11)
%     contourf(X, Y, NPPombsedgefrac,100,'LineStyle','none')
%     title('\fontsize{14}ombrotrophic sedge');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 colorscale])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
%     
%     subplot(4,3,12)
%     contourf(X, Y, NPPombshrubfrac,100,'LineStyle','none')
%     title('\fontsize{14}ombrotrohpic shrub');
%     xlim([-0.25 1])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 colorscale])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
% 
% end 
% 
% % plot or not
% 
% % plot_flag = 0;
% 
% if (plot_flag > 0)
% 
%   % ---------  
%   % figure(26)
%   % ---------  
% 
%     colorscale = 0.7;
%     figure(26)
%     subplot(2,1,1)
%     contourf(X, Y, NPPombsedgefrac + NPPminsedgefrac,100,'LineStyle','none')
%     title('\fontsize{14}total sedge NPPfrac');
% %    xlim([-0.25 1])
%     xlim([-0.1 0.5])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 colorscale])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
%     
%      subplot(2,1,2)
%     contourf(X, Y, squeeze((NPP(3,:,:)+NPP(11,:,:)))/npptotalmax,100,'LineStyle','none')
%     title('\fontsize{14}total sedge NPPfrac');
% %    xlim([-0.25 1])
%     xlim([-0.1 0.5])
%     ylim([0 8])
%     zlim([0 1])
%     colorbar
%     caxis([0 colorscale])
%     xlabel('\fontsize{14}water table depth [m]')
%     ylabel('\fontsize{14}peat depth [m]')
%     zlabel('\fontsize{14}relative NPP')
%     h1c=gca; 
%     set(h1c,'FontSize',14)
%     
   
% end

return