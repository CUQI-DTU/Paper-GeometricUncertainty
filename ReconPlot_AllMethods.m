%This script produces the plot of the reconstructions for the different
%methods (figure 6)

%Name of the figure
filename = 'my_figure_name';

figure('Position',[100 100 1000 1000])
cmin = 0; cmax = 0.13;

x_recon = zeros(2000^2,12);
cor_vec = zeros(1,12);
formatspec = '%0.2f';

%Load the results sequentially

%Load the high-dose results
load('my_highdose_results.mat');
x_recon(:,1) = res.X_MAP_ini; x_recon(:,2) = res.x_mean; x_recon(:,3) = res.X_MAP_com; x_recon(:,4) = res.X_MAP_xcorr;
cor_vec(1) = 0; cor_vec(2) = mean(res.cor_samps(end-1000:end)); cor_vec(3) = res.dx_com*res.setup.sy_true/(res.setup.sy_true + res.setup.dy_true); cor_vec(4) = res.dx_xcorr*res.setup.sy_true/(res.setup.sy_true + res.setup.dy_true);

%Load the low-dose results
load('my_lowdose_results.mat')
x_recon(:,5) = res.X_MAP_ini; x_recon(:,6) = res.x_mean; x_recon(:,7) = res.X_MAP_com; x_recon(:,8) = res.X_MAP_xcorr;
cor_vec(5) = 0; cor_vec(6) = mean(res.cor_samps(end-1000:end)); cor_vec(7) = res.dx_com*res.setup.sy_true/(res.setup.sy_true + res.setup.dy_true); cor_vec(8) = res.dx_xcorr*res.setup.sy_true/(res.setup.sy_true + res.setup.dy_true);

%Load the short-scan results
load('my_shortscan_results.mat');
x_recon(:,9) = res.X_MAP_ini; x_recon(:,10) = res.x_mean; x_recon(:,11) = res.X_MAP_com; x_recon(:,12) = res.X_MAP_xcorr;
cor_vec(9) = 0; cor_vec(10) = mean(res.cor_samps(end-1000:end)); cor_vec(11) = res.dx_com*res.setup.sy_true/(res.setup.sy_true + res.setup.dy_true); cor_vec(12) = res.dx_xcorr*res.setup.sy_true/(res.setup.sy_true + res.setup.dy_true);

%1200 and 400 are the default currently
region_x = 1100;
region_y = 500;

%cor_vec = dx_vec;
%ystart = [0.72, 0.51, 0.30, 0.09];
%ystart = [0.72,0.50,0.28,0.06] + 0.02;
ystart = [0.72,0.50,0.28,0.06];
xstart = [0.05 0.26 0.47 0.68];
%xstart = [0.05 0.26 0.57 0.68];
cor_vec = cor_vec/(17.4*10^(-3));

for i=1:3
    %ax1 = axes('Position',[xstart(5-i) ystart(1) 0.2 0.2]);
    %ax2 = axes('Position',[xstart(5-i) ystart(2) 0.2 0.2]);
    %ax3 = axes('Position',[xstart(5-i) ystart(3) 0.2 0.2]);
    
    ax1 = axes('Position',[xstart(i) ystart(1) 0.2 0.2]);
    ax2 = axes('Position',[xstart(i) ystart(2) 0.2 0.2]);
    ax3 = axes('Position',[xstart(i) ystart(3) 0.2 0.2]);
    ax4 = axes('Position',[xstart(i) ystart(4) 0.2 0.2]);
    im = reshape(x_recon(:,i*4),2000,2000);
    squarezoom_v2(im,4,250,region_x,region_y,ax4,[cmin,cmax]);
    %imagesc(ax1,im),caxis(ax1,[cmin cmax])
    title(ax4,['$c$ = ' num2str(cor_vec(i*4),formatspec)],'interpreter','latex')
    if i==1
        ylabel(ax4,'XCORR','interpreter','latex')
    end
    %xlabel(ax1,'energy bin')
    %ylabel(ax1,'detector $\times$ measurement')

    im = reshape(x_recon(:,i*4-1),2000,2000);
    squarezoom_v2(im,4,250,region_x,region_y,ax3,[cmin,cmax])
    title(ax3,['$c$ = ' num2str(cor_vec(i*4-1),formatspec)],'interpreter','latex')
    if i==1
        ylabel(ax3,'COM','interpreter','latex')
    end
    %xlabel(ax2,'energy bin')
    %ylabel(ax2,'detector $\times$ measurement')

    im = reshape(x_recon(:,i*4-2),2000,2000);
    squarezoom_v2(im,4,250,region_x,region_y,ax2,[cmin,cmax])
    %title(ax3,'(a) image 3')
    %xlabel(ax3,'energy bin')
    %ylabel(ax3,'detector $\times$ measurement')
    title(ax2,['$c$ = ' num2str(cor_vec(i*4-2),formatspec)],'interpreter','latex')
    if i==1
        ylabel(ax2,'MCMC','interpreter','latex')
    end
    
    
    im = reshape(x_recon(:,i*4-3),2000,2000);
    squarezoom_v2(im,4,250,region_x,region_y,ax1,[cmin,cmax])
    if i == 1
        title(ax1,{'High-dose','$c$ = 0'},'interpreter','latex')
        ylabel(ax1,'Initial','interpreter','latex')
    elseif i==2
        title(ax1, {'Low-dose','$c$ = 0'},'interpreter','latex')
    elseif i==3
        title(ax1, {'Short-scan','$c$ = 0'},'interpreter','latex')
    end
    
    set(ax1, 'XTick', []);
    set(ax1, 'YTick', []);
    set(ax2, 'XTick', []);
    set(ax2, 'YTick', []);
    set(ax3, 'XTick', []);
    set(ax3, 'YTick', []);
    set(ax4, 'XTick', []);
    set(ax4, 'YTick', []);
end
colorbar(ax4,'location','southoutside','Position',...
      [0.05 0.03 0.62 0.02]);

%Save figure
saveas(gcf,filename,'epsc')
