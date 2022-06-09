function MCMC_Plot_Parallel(res,N_burn,N_trim,int_recon,int_cred)

%Unpack variables
x_samps = res.x_samps;
delta_samps = res.delta_samps;
lambda_samps = res.lambda_samps;

model_param_count = res.setup.DETECTOR_X;

figure(1)
h1 = tight_subplot(model_param_count+2,3,[0.09,0.04],[0.07,0.05]);
figure(2)
h2 = tight_subplot(model_param_count+2,4,[0.09,0.04],[0.07,0.05]);

axes(h1(1))
plot(1:length(lambda_samps), lambda_samps)
title('\lambda - chain')
axes(h1(2))
histogram(lambda_samps)
title('\lambda - histogram')
axes(h1(3))
autocorr(lambda_samps)
title('\lambda - ACF')
ylabel('')
axes(h1(4))
plot(1:length(delta_samps), delta_samps)
title('\delta - chain')
axes(h1(5))
histogram(delta_samps)
title('\delta - histogram')
axes(h1(6))
autocorr(delta_samps)
ylabel('')
title('\delta - ACF')

figure(2)
axes(h2(1))
plot(1:length(lambda_samps),lambda_samps)
title('\lambda - chain - full')
xlim([1,length(lambda_samps)])
ylim([min(lambda_samps)-1000,max(lambda_samps)+1000])
axes(h2(2))
plot(lambda_samps(N_burn+1:N_trim:end))
title('\lambda - chain')
axes(h2(3))
histogram(lambda_samps(N_burn+1:N_trim:end))
title('\lambda - histogram')
axes(h2(4))
autocorr(lambda_samps(N_burn+1:N_trim:end))
title('\lambda - ACF')
ylabel('')
axes(h2(5))
plot(1:length(delta_samps),delta_samps)
xlim([1,length(delta_samps)])
ylim([0,max(delta_samps)+10])
title('\delta - chain - full')
axes(h2(6))
plot(delta_samps(N_burn+1:N_trim:end))
title('\delta - chain')
axes(h2(7))
histogram(delta_samps(N_burn+1:N_trim:end))
title('\delta - histogram')
axes(h2(8))
autocorr(delta_samps(N_burn+1:N_trim:end))
ylabel('')
title('\delta - ACF')

mod_count = 0;

if res.setup.DETECTOR_X == 1
    dx_true = res.setup.dx_true;
    dx_samps = res.dx_samps;
    mod_count = mod_count + 1;
    
    axes(h1(6+mod_count*3-2))
    plot(1:length(dx_samps),dx_samps)
    title('x_d-chain')
    axis tight
    axes(h1(6+mod_count*3-1))
    histogram(dx_samps)
    title('x_d-histogram')
    axes(h1(6+mod_count*3))
    autocorr(dx_samps)
    ylabel('')
    title('x_d-ACF')
        
    axes(h2(8+mod_count*4-3))
    plot(dx_samps)
    title('d - chain - full')
    xlim([1,length(dx_samps)])
    if dx_true<min(dx_samps)
        ylim([dx_true,max(dx_samps)+0.1])
    elseif dx_true>max(dx_samps)
        ylim([min(dx_samps),dx_true+0.1])
    else
        ylim([min(dx_samps),max(dx_samps)+0.1])
    end
    yline(dx_true,'--r');
    axes(h2(8+mod_count*4-2))
    plot(1:(length(dx_samps)-N_burn),dx_samps(N_burn+1:end))
    title('d - chain')
    axes(h2(8+mod_count*4-1))
    histogram(dx_samps(N_burn+1:end))
    title('d - histogram')
    axes(h2(8+mod_count*4))
    autocorr(dx_samps(N_burn+1:end))
    ylabel('')
    title('d - ACF')
    
   % n_accept_dx = res.n_accept_dx;
   % acceptance_rate_dx = sum(n_accept_dx(N_burn+1:end))/((res.setup.N_samples-N_burn)*res.setup.N_metro);
    
   % fprintf('dx true value: %f\n',dx_true)
   % fprintf('dx initial value: %f\n',res.setup.dx0)
   % fprintf('dx MH step size: %f\n',res.setup.dx_sigma_proposal)
   % fprintf('dx acceptance rate: %f\n\n',acceptance_rate_dx)
end

N = res.setup.N;

if res.setup.nonneg == 1
    con = 'non-negative';
    fprintf('Non-negativity Constraints Enforced\n')
else
    con = '';
    fprintf('No non-negativity Constraints Enforced\n')
end

if res.setup.iid == 1
    reg_term = 'Tikh.';
    fprintf('Ordinary Tikhonov Regularization\n\n')
else
    reg_term = 'Gentikh.';
    fprintf('Generalized Tokhonov Regularization\n\n')
end

%Plot sample mean reconstruction and width of 95% pixelwise credibility
%Remove Burn in samples
x_samps(:,1:N_burn) = [];

x_mean = mean(x_samps(:,1:N_trim:end),2);
lower_quant = quantile(x_samps,0.025,2);
upper_quant = quantile(x_samps,0.975,2);

if nargin>2
    figure
    h3 = tight_subplot(1,2);
    axes(h3(1))
    imagesc(reshape(x_mean,N,N),int_recon) 
    axis image
    colorbar
    title('Sample Mean')
    axes(h3(2))
    imagesc(reshape(upper_quant-lower_quant,N,N),int_cred)
    axis image
    colorbar
    title('Width of 95% Credibility Interval')
    
    figure
    h4 = tight_subplot(2,2,[0.1,0.02],[0.05,0.05],[0.05,0.01]);
    axes(h4(1))
    imagesc(reshape(res.X_MAP_ini,N,N),int_recon)
    axis image
    colorbar
    title('MAP - Initial Geometry')
    axes(h4(2))
    imagesc(reshape(res.X_MAP_true,N,N),int_recon)
    axis image
    colorbar
    title('MAP - True Geometry')
    axes(h4(3))
    imagesc(reshape(x_mean,N,N),int_recon)
    axis image
    colorbar
    title('Sample Mean')
    axes(h4(4))
    imagesc(reshape(upper_quant-lower_quant,N,N),int_cred)
    axis image
    colorbar
    title('Width of 95% Credibility Interval')
    
    figure
    h5 = tight_subplot(1,3,[0.1,0.02],[0.05,0.05],[0.05,0.01]);
    axes(h5(1))
    imagesc(reshape(x_mean,N,N),int_recon)
    title('Sample Mean')
    axes(h5(2))
    imagesc(reshape(X_MAP_com,N,N),int_recon)
    title('COM')
    axes(h5(3))
    imagesc(reshape(X_MAP_xcorr,N,N),int_recon)
    title('XCORR')

    figure
    imagesc(reshape(res.X_MAP_ini,N,N),int_recon), axis image, colorbar
    title('Initial MAP reconstruction')

    figure
    imagesc(reshape(res.setup.x0,N,N),int_recon), axis image, colorbar
    title('Initial MAP Reconstruction - wrong reg')

    figure
    imagesc(reshape(x_mean,N,N),int_recon), axis image, colorbar
    title('Sample Mean')

    figure
    imagesc(reshape(upper_quant-lower_quant,N,N),int_cred), axis image, colorbar
    title('Width of 95% Credibility Interval')
else
    figure
    h3 = tight_subplot(1,2);
    axes(h3(1))
    imagesc(reshape(x_mean,N,N))
    axis image
    colorbar
    title('Sample Mean')
    axes(h3(2))
    imagesc(reshape(upper_quant-lower_quant,N,N))
    axis image
    colorbar
    title('Width of 95% Credibility Interval')
    
    figure
    h4 = tight_subplot(2,2,[0.1,0.05],[0.05,0.05],[0.05,0.01]);
    axes(h4(1))
    imagesc(reshape(res.setup.x0,N,N))
    axis image
    colorbar
    title('MAP - Initial Geometry')
    axes(h4(2))
    imagesc(reshape(res.X_MAP_mean,N,N))
    axis image
    colorbar
    title('MAP - Center of Mass Alignment')
    axes(h4(3))
    imagesc(reshape(x_mean,N,N))
    axis image
    colorbar
    title('Sample Mean')
    axes(h4(4))
    imagesc(reshape(upper_quant-lower_quant,N,N))
    axis image
    colorbar
    title('Width of 95% Credibility Interval')

    figure
    imagesc(reshape(res.X_MAP_mean,N,N)), axis image, colorbar
    title('Initial MAP reconstruction')

    figure
    imagesc(reshape(res.setup.x0,N,N)), axis image, colorbar
    title('Initial MAP Reconstruction - wrong reg')

    figure
    imagesc(reshape(x_mean,N,N)), axis image, colorbar
    title('Sample Mean')

    figure
    imagesc(reshape(upper_quant-lower_quant,N,N)), axis image, colorbar
    title('Width of 95% Credibility Interval')
end

%Plot images of upper and lower quantiles
figure
subplot(1,2,1)
imagesc(reshape(lower_quant,N,N),[min(min(lower_quant)), max(max(upper_quant))]), axis image, colorbar
title('Lower Quantile')
subplot(1,2,2)
imagesc(reshape(upper_quant,N,N),[min(min(lower_quant)), max(max(upper_quant))]), axis image, colorbar
title('Upper Quantile')

figure
imagesc(reshape(res.setup.b,length(res.setup.theta),res.setup.p)), axis image, colorbar
title('Noisy Sinogram')

%Do geweke test
[z_lambda,p_lambda] = geweke(lambda_samps(N_burn:end)');
[z_delta,p_delta] = geweke(delta_samps(N_burn:end)');
disp(['p-value for lambda chain: ' num2str(p_lambda)])
disp(['p-value for delta chain: ' num2str(p_delta)])

if res.setup.DETECTOR_X == 1
    [z_dx,p_dx] = geweke(dx_samps(N_burn:end)');
    disp(['p-value for dx-chain: ' num2str(p_dx)])
end

end