%This script produces the reconstruction plot (figure 2)
%%
figure('Position',[100 100 1000 1000])
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultColorbarTickLabelInterpreter','latex');

cmin = 0; cmax = 0.13;

%Load results
load('COR_BeadsRecon_MAP.mat');

%cor_vec = dx_vec;
%ystart = [0.72, 0.51, 0.30, 0.09];
%ystart = [0.72,0.50,0.28,0.06] + 0.02;
ystart = [0.72,0.50,0.28,0.06];
xstart = [0.05 0.26 0.47 0.68];
%xstart = [0.05 0.26 0.57 0.68];

%Convert to pixel units
cor_vec = cor_vec/(17.4*10^(-3));

region_x = 1100;
region_y = 500;

for i=1:2
    %ax1 = axes('Position',[xstart(5-i) ystart(1) 0.2 0.2]);
    %ax2 = axes('Position',[xstart(5-i) ystart(2) 0.2 0.2]);
    %ax3 = axes('Position',[xstart(5-i) ystart(3) 0.2 0.2]);
    
    ax1 = axes('Position',[xstart(1) ystart(i) 0.2 0.2]);
    ax2 = axes('Position',[xstart(2) ystart(i) 0.2 0.2]);
    ax3 = axes('Position',[xstart(3) ystart(i) 0.2 0.2]);
    im = reshape(x_recon(:,i*3),2000,2000);
    squarezoom_v2(im,4,250,region_x,region_y,ax3,[cmin,cmax]);
    %imagesc(ax1,im),caxis(ax1,[cmin cmax])
    title(ax3,['c = ' num2str(cor_vec(i*3))],'interpreter','latex')
    %xlabel(ax1,'energy bin')
    %ylabel(ax1,'detector $\times$ measurement')

    im = reshape(x_recon(:,i*3-1),2000,2000);
    squarezoom_v2(im,4,250,region_x,region_y,ax2,[cmin,cmax])
    title(ax2,['c = ' num2str(cor_vec(i*3-1))],'interpreter','latex')
    %xlabel(ax2,'energy bin')
    %ylabel(ax2,'detector $\times$ measurement')

    im = reshape(x_recon(:,i*3-2),2000,2000);
    squarezoom_v2(im,4,250,region_x,region_y,ax1,[cmin,cmax])
    %title(ax3,'(a) image 3')
    %xlabel(ax3,'energy bin')
    %ylabel(ax3,'detector $\times$ measurement')
    title(ax1,['c = ' num2str(cor_vec(i*3-2))],'interpreter','latex')
    
    set(ax1, 'XTick', []);
    set(ax1, 'YTick', []);
    set(ax2, 'XTick', []);
    set(ax2, 'YTick', []);
    set(ax3, 'XTick', []);
    set(ax3, 'YTick', []);
end
colorbar(ax3,'location','southoutside','Position',...
      [0.05 0.47 0.62 0.02]);
  
saveas(gcf,'ReconplotCOR','epsc')
