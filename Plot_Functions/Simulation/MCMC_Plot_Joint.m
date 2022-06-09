function MCMC_Plot_Joint(res,N_burn,N_trim,int_recon,int_cred)

%Unpack variables
x_samps = res.x_samps;
delta_samps = res.delta_samps;
lambda_samps = res.lambda_samps;

model_param_count = res.setup.SOURCE_X + res.setup.SOURCE_Y + res.setup.DETECTOR_X + res.setup.DETECTOR_Y;

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

if res.setup.SOURCE_X == 1
    sx_true = res.setup.sx_true;
    sx_samps = res.sx_samps;
    mod_count = mod_count + 1;
    axes(h1(6+mod_count*3-2))
    plot(1:length(sx_samps),sx_samps)
    title('x_s-chain')
    axes(h1(6+mod_count*3-1))
    histogram(sx_samps)
    title('x_s-histogram')
    axes(h1(6+mod_count*3))
    autocorr(sx_samps)
    ylabel('')
    title('x_s - ACF')
        
    axes(h2(8+mod_count*4-3))
    plot(sx_samps)
    title('x_s - chain - full')
    xlim([1,length(sx_samps)])
    if sx_true<min(sx_samps)
        ylim([sx_true,max(sx_samps)+0.1])
    elseif sx_true>max(sx_samps)
        ylim([min(sx_samps),sx_true+0.1])
    else
        ylim([min(sx_samps),max(sx_samps)+0.1])
    end
    yline(sx_true,'--r');
    axes(h2(8+mod_count*4-2))
    plot(sx_samps(N_burn+1:N_trim:end))
    title('x_s - chain')
    axes(h2(8+mod_count*4-1))
    histogram(sx_samps(N_burn+1:N_trim:end))
    title('x_s - histogram')
    axes(h2(8+mod_count*4))
    autocorr(sx_samps(N_burn+1:N_trim:end))
    ylabel('')
    title('x_s - ACF')
    
    fprintf('sx true value: %f\n',sx_true)
    fprintf('sx initial value: %f\n',res.setup.sx0)
    fprintf('sx MH step size: %f\n\n',res.setup.sx_sigma_proposal)
end

if res.setup.SOURCE_Y == 1
    sy_true = res.setup.sy_true;
    sy_samps = res.sy_samps;
    mod_count = mod_count + 1;
    
    axes(h1(6+mod_count*3-2))
    plot(1:length(sy_samps),sy_samps)
    title('y_s - chain - full')
    yline(sy_true);
    axes(h1(6+mod_count*3-1))
    histogram(sy_samps)
    title('y_s - histogram')
    axes(h1(6+mod_count*3))
    autocorr(sy_samps)
    ylabel('')
    title('y_s - ACF')
        
    axes(h2(8+mod_count*4-3))
    plot(sy_samps)
    title('y_s - chain - full')
    xlim([1,length(sy_samps)])
    if sy_true<min(sy_samps)
        ylim([sy_true,max(sy_samps)+0.5])
    elseif sy_true>max(sy_samps)
        ylim([min(sy_samps),sy_true+0.5])
    else
        ylim([min(sy_samps),max(sy_samps)+0.5])
    end
    yline(sy_true,'--r');
    axes(h2(8+mod_count*4-2))
    plot(sy_samps(N_burn+1:N_trim:end))
    title('y_s - chain')
    axes(h2(8+mod_count*4-1))
    histogram(sy_samps(N_burn+1:N_trim:end))
    title('y_s - histogram')
    axes(h2(8+mod_count*4))
    autocorr(sy_samps(N_burn+1:N_trim:end))
    ylabel('')
    title('y_s - ACF')
    
    
    fprintf('sy true value: %f\n',sy_true)
    fprintf('sy initial value: %f\n',res.setup.sy0)
    fprintf('sy MH step size: %f\n\n',res.setup.sy_sigma_proposal)
end

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
    title('x_d - chain - full')
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
    plot(dx_samps(N_burn+1:N_trim:end))
    title('x_d - chain')
    axes(h2(8+mod_count*4-1))
    histogram(dx_samps(N_burn+1:N_trim:end))
    title('x_d - histogram')
    axes(h2(8+mod_count*4))
    autocorr(dx_samps(N_burn+1:N_trim:end))
    ylabel('')
    title('x_d - ACF')
    
    
    fprintf('dx true value: %f\n',dx_true)
    fprintf('dx initial value: %f\n',res.setup.dx0)
    fprintf('dx MH step size: %f\n\n',res.setup.dx_sigma_proposal)
end

if res.setup.DETECTOR_Y == 1
    dy_true = res.setup.dy_true;
    dy_samps = res.dy_samps;
    mod_count = mod_count + 1;
    
    axes(h1(6+mod_count*3-2))
    plot(1:length(dy_samps),dy_samps)
    title('y_d - chain')
    axes(h1(6+mod_count*3-1))
    histogram(dy_samps)
    title('y_d - histogram')
    axes(h1(6+mod_count*3))
    autocorr(dy_samps)
    ylabel('')
    title('y_d - ACF')
        
    axes(h2(8+mod_count*4-3))
    plot(dy_samps)
    title('y_d - chain - full')
    xlim([1,length(dy_samps)])
    if dy_true<min(dy_samps)
        ylim([dy_true,max(dy_samps)+0.5])
    elseif dy_true>max(dy_samps)
        ylim([min(dy_samps),dy_true+0.5])
    else
        ylim([min(dy_samps),max(dy_samps)+0.5])
    end
    yline(dy_true,'--r');
    axes(h2(8+mod_count*4-2))
    plot(dy_samps(N_burn+1:N_trim:end))
    title('y_d - chain')
    axes(h2(8+mod_count*4-1))
    histogram(dy_samps(N_burn+1:N_trim:end))
    title('y_d - histogram')
    axes(h2(8+mod_count*4))
    autocorr(dy_samps(N_burn+1:N_trim:end))
    ylabel('')
    title('y_d - ACF')
    
    fprintf('dy true value: %f\n',dy_true)
    fprintf('dy initial value: %f\n',res.setup.dy0)
    fprintf('dy MH step size: %f\n\n',res.setup.dy_sigma_proposal)
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

accept_rate = sum(res.accept_vec(N_burn:end))/((res.setup.N_samples-N_burn)*res.setup.N_metro);
fprintf('MH acceptance rate: %f\n',accept_rate)
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

if res.setup.SOURCE_X == 1
    [z_sx,p_sx] = geweke(sx_samps(N_burn:end)');
    disp(['p-value for sx-chain: ' num2str(p_sx)])
end

if res.setup.SOURCE_Y == 1
    [z_sy,p_sy] = geweke(sy_samps(N_burn:end)');
    disp(['p-value for sy-chain: ' num2str(p_sy)])
end

if res.setup.DETECTOR_X == 1
    [z_dx,p_dx] = geweke(dx_samps(N_burn:end)');
    disp(['p-value for dx-chain: ' num2str(p_dx)])
end

if res.setup.DETECTOR_Y == 1
    [z_dy,p_dy] = geweke(dy_samps(N_burn:end)');
    disp(['p-value for dy-chain: ' num2str(p_dy)])
end
end