%This script produces the chain plots (figure 3,4,5) in the article.

%Load results file produced by 'Simple_MCMC_driver.m' script
load('my_results.mat')

%filename for produced plot
filename = 'my_chain_plot';

%1000 600 first
figure('Position',[100 100 1000 600])

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultColorbarTickLabelInterpreter','latex');
%formatspec = '%0.3f';

%%

%ystart = [0.72, 0.51, 0.30, 0.09];
%ystart = [0.72,0.50,0.28,0.06] + 0.02;
ystart = [0.68,0.36,0.04];
xstart = [0.04 0.295 0.55 0.795]-0.005;
%xstart = [0.05 0.26 0.57 0.68];

%set(0,'defaultTextInterpreter','latex');
N_burn = res.setup.Nburn;
N_trim = 1;
for i=1:3
    %ax1 = axes('Position',[xstart(5-i) ystart(1) 0.2 0.2]);
    %ax2 = axes('Position',[xstart(5-i) ystart(2) 0.2 0.2]);
    %ax3 = axes('Position',[xstart(5-i) ystart(3) 0.2 0.2]);
    
    ax1 = axes('Position',[xstart(1) ystart(i) 0.2 0.25]);
    ax2 = axes('Position',[xstart(2) ystart(i) 0.2 0.25]);
    ax3 = axes('Position',[xstart(3) ystart(i) 0.2 0.25]);
    ax4 = axes('Position',[xstart(4) ystart(i) 0.2 0.25]);
    
    if i==1
        chain = res.lambda_samps;
        colorstring = '#D95319';
        colortriplet = [0.8500 0.3250 0.0980];
        name = '$\lambda$';
        offset = 1000;
    elseif i==2
        chain = res.delta_samps;
        colorstring = '#0072BD';
        colortriplet = [0 0.4470 0.7410];
        name = '$\delta$';
        offset = 50;
    elseif i==3
        chain = res.cor_samps/(17.4*10^(-3));
        colorstring = '#77AC30';
        colortriplet = [0.4660 0.6740 0.1880];
        name = '$c$';
        offset = 1;
    end
    
    q_lower = quantile(chain(N_burn:end),0.025);
    q_upper = quantile(chain(N_burn:end),0.975);
    
    plot(ax1,1:length(chain),chain,'Color',colorstring)
    title(ax1,[name ' - chain - full'])
    xlim(ax1,[1,length(chain)])
    ylim(ax1,[min(chain)-offset,max(chain)+offset])
    if i==3
        ax1.YTick = [0,3,6,9,12];
    end
    plot(ax2,(N_burn+1):length(chain),chain(N_burn+1:N_trim:end),'Color',colorstring)
    yline(ax2,q_lower,'--k')
    yline(ax2,q_upper,'--k')
    title(ax2,[name ' - chain'])
    if i==3
         %ax2.YAxis.TickLabelFormat = '%.4f';
         %ax2.XLim = [12.240,12.248];
    end
    histogram(ax3,chain(N_burn+1:N_trim:end),'Normalization','probability','Facecolor',colorstring)
    xline(ax3,q_lower,'--k')
    xline(ax3,q_upper,'--k')
    title(ax3,[name ' - histogram'])
    [~,~,~,acf] = autocorr(chain(N_burn+1:N_trim:end));
    acf(1).Color = colortriplet;
    acf(2).Color = [0 0 0];
    acf(3).Color = [0 0 0];
    acf(2).LineStyle = '--';
    acf(3).LineStyle = '--';
    title(ax4,[name ' - ACF'])
    xlabel(ax4,[])
    ylabel(ax4,[])
end
saveas(gcf,filename,'epsc')
