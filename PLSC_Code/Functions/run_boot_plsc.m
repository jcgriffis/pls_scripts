function [plsc] = run_boot_plsc(plsc, cfg, sig_comps, n_obs)
% This function computes runs bootstraps on the PLSC between a predictor matrix and a response
% matrix.

% Inputs:
% plsc structure created by run_plsc and augmented by run_b_plsc
% cfg structure created prior to running run_plsc
% sig_comps: indices of significant components based on p-values obtained from run_b_plsc
% n_obs: number of observations

% Outputs: 
% Augmented plsc structure with bootstrap results

% Joseph Griffis 2018, Washington University in St Louis

disp('Running bootstraps...');
plsc.boot_r = zeros(cfg.n_boot, length(sig_comps));
for j = 1:cfg.n_boot
    disp(['Boostrap number ' num2str(j) ' of ' num2str(cfg.n_boot)]);
    [X_boot, idx] = datasample(cfg.X, n_obs, 1, 'Replace', true); % resample X
    Y_boot = cfg.Y(idx,:);
    
    % center and re-scale X and Y and remove NaNs
    if cfg.zscore == 0
        nX_b = X_boot - (repmat(mean(X_boot),[n_obs,1]));
        nY_b = Y_boot - (repmat(mean(Y_boot),[n_obs,1]));
    elseif cfg.zscore == 1
        disp('Z-score flag set to 1: Standardizing variables...');
        % rescale variables
        nX_b = X_boot - (repmat(mean(X_boot),[n_obs,1]))./(repmat(std(X_boot),[n_obs,1]));
        nY_b = Y_boot - (repmat(mean(Y_boot),[n_obs,1]))./(repmat(std(Y_boot),[n_obs,1]));
        nX_b(isnan(nX_b))=0; nX_b(isinf(nX_b))=0;
        nY_b(isnan(nY_b))=0; nY_b(isinf(nY_b))=0;
    end
    
    % PLSC on resampled data
    R_boot = nY_b'*nX_b; % get covariance of X and Y
    [U_b, S_b, V_b] = paq(R_boot, length(diag(plsc.S))); % compute SVD
    if length(S_b) < length(diag(plsc.S))
        U_b = [U_b, zeros(size(U_b,1), length(diag(plsc.S))-length(S_b))];
        V_b = [V_b, zeros(size(V_b,1), length(diag(plsc.S))-length(S_b))];
        S_b = [S_b; zeros(length(diag(plsc.S))-length(S_b),1)];
    end
    
    S_b = diag(S_b); % get singular values
    
    [N, O, P] = paq(plsc.V'*V_b); % procruste analysis to match LVs
    
    Q = N*P'; % get transform
    V_hat_b = V_b*S_b*Q; % apply transform to V
    U_hat_b = U_b*S_b*Q; % apply transform to U
     
    % Get bootstrapped statistics for significant components
    boot_U(j,:,1:length(sig_comps)) = U_hat_b(:,sig_comps); % statistics for U
    boot_V(j,:,1:length(sig_comps)) = V_hat_b(:,sig_comps); % statistics for V
    
    % get correlations
    Lx_b = nX_b*V_hat_b; % get latent variable scores
    Ly_b = nY_b*U_hat_b; % get latent variable scores
    for i = 1:length(sig_comps)
        plsc.boot_r(j,i) = corr(Lx_b(:,i), Ly_b(:,i), 'rows', 'pairwise');
    end
    
end

disp('Bootstrap analyses complete');
plsc.bs_ratio_Y = [];
plsc.bs_ratio_X = [];
for i = 1:length(sig_comps)
    se_u = std(boot_U(:,:,i),1)'; % get SD 
    plsc.bs_ratio_Y(:,i) = plsc.Fy(:,sig_comps(i))./se_u; % get ratio    
    se_v = std(boot_V(:,:,i),1)'; % get SD 
    plsc.bs_ratio_X(:,i) = plsc.Fx(:,sig_comps(i))./se_v; % get ratio 
end
end