function Model_Discrepency_Plot(res,int_recon)
%This function plots parameter estimates and reconstructions from the
%Model discrepency algorithm

out = res.out;
setup = res.setup;

reg_term = setup.reg_term;
mod_param_count = setup.DETECTOR_X + setup.DETECTOR_Y + setup.SOURCE_X + setup.SOURCE_Y;

if setup.nonneg == 1
    string = 'nonneg';
else
    string = '';
end

x_final_recon = out.x_final_recon;
x_initial_recon = out.x_initial_recon;
x_final_MAP = out.x_final_MAP;

mod_count = 0;
    
N = sqrt(length(x_final_recon));
figure(1)
if setup.DETECTOR_X == 1
    dx0 = setup.dx0;
    dx_true = setup.dx_true;
    dx_samps = out.dx_samps;
    mod_count = mod_count + 1; 
    if mod_param_count == 4
        figure(1)
        subplot(2,2,mod_count)
        plot(0:length(dx_samps),[dx0,dx_samps],'x-')
        yline(dx_true,'--r')
        ylim([min([dx_samps,dx_true])-0.1,max([dx_samps,dx_true])+0.1])
        title('dx - Parameter Estimates')
    else
        figure(1)
        subplot(mod_param_count,1,mod_count)
        plot(0:length(dx_samps),[dx0,dx_samps],'x-')
        yline(dx_true,'--r')
        ylim([min([dx_samps,dx_true])-0.1,max([dx_samps,dx_true])+0.1])
        title('dx - Parameter Estimates')
    end
end
        
if setup.DETECTOR_Y == 1
    dy0 = setup.dy0;
    dy_true = setup.dy_true;
    dy_samps = out.dy_samps;
    mod_count = mod_count + 1;
    if mod_param_count == 4
        figure(1)
        subplot(2,2,mod_count)
        plot(0:length(dy_samps),[dy0,dy_samps],'x-')
        yline(dy_true,'--r')
        ylim([min([dy_samps,dy_true])-1,max([dy_samps,dy_true])+1])
        title('dy - Parameter Estimates')
    else
        figure(1)
        subplot(mod_param_count,1,mod_count)
        plot(0:length(dy_samps),[dy0,dy_samps],'x-')
        yline(dy_true,'--r')
        ylim([min([dy_samps,dy_true])-1,max([dy_samps,dy_true])+1])
        title('dy - Parameter Estimates')
    end
end
if setup.SOURCE_X == 1
    sx0 = setup.sx0;
    sx_true = setup.sx_true;
    sx_samps = out.sx_samps;
    mod_count = mod_count + 1;
    if mod_param_count == 4
        figure(1)
        subplot(2,2,mod_count)
        plot(0:length(sx_samps),[sx0,sx_samps],'x-')
        yline(sx_true,'--r')
        ylim([min([sx_samps,sx_true])-0.2,max([sx_samps,sx_true])+0.2])
        title('sx - Parameter Estimates')
    else
        figure(1)
        subplot(mod_param_count,1,mod_count)
        plot(0:length(sx_samps),[sx0,sx_samps],'x-')
        yline(sx_true,'--r')
        ylim([min([sx_samps,sx_true])-0.2,max([sx_samps,sx_true])+0.2])
        title('sx - Parameter Estimates')
    end
end
if setup.SOURCE_Y == 1
    sy0 = setup.sy0;
    sy_true = setup.sy_true;
    sy_samps = out.sy_samps;
    mod_count = mod_count + 1;
    if mod_param_count == 4
        figure(1)
        subplot(2,2,mod_count)
        plot(0:length(sy_samps),[sy0,sy_samps],'x-')
        yline(sy_true,'--r')
        ylim([min([sy_samps,sy_true])-1,max([sy_samps,sy_true])+1])
        title('sy - Parameter Estimates')
    else
        figure(1)
        subplot(mod_param_count,1,mod_count)
        plot(0:length(sy_samps),[sy0,sy_samps],'x-')
        yline(sy_true,'--r')
        ylim([min([sy_samps,sy_true])-1,max([sy_samps,sy_true])+1])
        title('sy - Parameter Estimates')
    end
end 
    
%Plot the initial reconstruction, the final outer reconstruction and the
%reconstruction with the true model parameters
if nargin>1
    figure
    subplot(1,3,1)
    imagesc(reshape(x_initial_recon,N,N),int_recon), axis image, colorbar
    title('MAP - Initial Parameters')
    subplot(1,3,2)
    imagesc(reshape(x_final_MAP,N,N),int_recon), axis image, colorbar
    title('MAP - Final Parameter Estimates')
    subplot(1,3,3)
    imagesc(reshape(x_true_model,N,N),int_recon), axis image, colorbar
    title('MAP - True Parameters')

    figure
    imagesc(reshape(x_initial_recon,N,N),int_recon), axis image, colorbar
    title('MAP - Initial Parameters')

    figure
    imagesc(reshape(x_final_MAP,N,N),int_recon), axis image, colorbar
    title('MAP - Final Parameter Estimates')

    figure
    imagesc(reshape(x_final_recon,N,N),int_recon), axis image, colorbar
    title('Final Outer Reconstruction')
else
    figure
    subplot(1,3,1)
    imagesc(reshape(x_initial_recon,N,N)), axis image, colorbar
    title('MAP - Initial Parameters')
    subplot(1,3,2)
    imagesc(reshape(x_final_MAP,N,N)), axis image, colorbar
    title('MAP - Final Parameter Estimates')
    subplot(1,3,3)
    imagesc(reshape(x_true_model,N,N)), axis image, colorbar
    title('MAP - True Parameters')

    figure
    imagesc(reshape(x_initial_recon,N,N)), axis image, colorbar
    title('MAP - Initial Parameters')

    figure
    imagesc(reshape(x_final_MAP,N,N)), axis image, colorbar
    title('MAP - Final Parameter Estimates')

    figure
    imagesc(reshape(x_final_recon,N,N)), axis image, colorbar
    title('Final Outer Reconstruction')
end