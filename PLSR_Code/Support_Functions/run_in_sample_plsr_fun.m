function [plsr_results] = run_in_sample_plsr_fun(cfg)
% This function runs in-sample PLSR analyses to predict a response variable (e.g. behavior) from an imaging predictor matrix
% and uses leave-one-out cross-validation to determine the optimal number of predictor components. 
% Bootstrapping is used to obtain CIs on model fit and betas (i.e. see Griffis et al., 2019, Cell Reports).
% Joseph Griffis, 2020

%%%%%%% Inputs %%%%%%%%
% X - an observation-by-predictor (e.g. patient-by-voxel, patient-by-edge), predictor matrix (rows are observations, columns are variables).
% Y - a column vector of responses to be used as the dependent variable.
% trim_flag - a flag for trimming predictor matrix to remove columns with insufficient observations (e.g. voxels that are only lesioned in 1 patient).
% freq_thresh - frequency threshold (i.e. how many patients need non-zero data for a given predictor column to keep it; only used if trim_flag == 1)
% alpha_thresh - desired alpha threshold
% fwe_flag - flag for FWE correction
% n_boot - number of bootstraps to run for CI calculation (e.g. 1000)

%%%%%% Outputs %%%%%%%%
% plsr_results - structure containing relevant results and parameters (see end of function)


% trim predictor data if flagged
if cfg.trim_flag == 1
    disp(['Trim flag set to 1, trimming data to retain columns with at least ' num2str(cfg.freq_thresh) ' observations greater than ' num2str(cfg.min_obs)]);
    [X_trim, X_inc_ind] = get_has_dmg(cfg.X(cfg.include,:), cfg.min_obs, cfg.freq_thresh); % get trimmed X_trim matrix and indices for retained cells
else
    disp(['Trim flag set to 0, using entire predictor matrix']);
    X_trim = cfg.X(cfg.include,:);
    X_inc_ind = 1:length(X_trim);
end
X_trim(isnan(X_trim)==1)=0; % remove any NANs from X_trim
X_trim(isinf(X_trim)==1)=0; % remove any Infs from X_trim
Y = cfg.Y(cfg.include); % final response variable

% run LOOCV optimmization for number of PLS components
disp('Running LOOCV optimization for PLS components');
[opt_k] = run_loocv_plsr(X_trim, Y); % loocv prediction to identify opt_k for fixed effects and get out-of-sample prediction
disp(['The optimal number of PLS components is ' num2str(opt_k) '. Fitting full model with optimal components.']);
[~,~,XS,~,beta] = plsregress(X_trim, Y, opt_k); % fit fixed effects model with opt_k
yhat = [ones(length(Y), 1) X_trim]*beta; % get fitted Y
numerator = nansum((Y-yhat).^2); % get squared error
denominator = nansum((Y-nanmean(Y)).^2);
r2 = 1-(numerator/denominator);

disp(['The full model explains ' num2str(100*r2) ' % of the variance in the response.']);

% Get bootstrapped confidence intervals for disconnection models
if cfg.fwe_flag == 1
    alpha_thresh = cfg.alpha_thresh./size(X_trim,2); % family-wise error corrected
    disp(['FWE flag set to 1: running boostrapped models to obtain ' num2str((1-cfg.alpha_thresh).*100) '% FWE-corrected confidence intervals']); 
else
    disp(['FWE flag set to 0: running boostrapped models to obtain ' num2str((1-cfg.alpha_thresh).*100) '% uncorrected confidence intervals']); 
end
ci_type = cfg.ci_type;
rng_seed = rng; % rng seed
boot_plsr_fun = @(IV,DV,k) run_boot_plsr(IV,DV,k); % returns R^2, AIC, and betas as a vector
options = statset('UseParallel', 'always'); % parallel processing
[ci] = bootci(cfg.n_boot, {boot_plsr_fun, X_trim, Y, opt_k}, 'alpha', alpha_thresh, 'Options', options, 'type', ci_type); % run bootstrapped PLSRs to get R^2, AIC, and betas

% threshold betas to retain non-zero crossing CIs
disp(['Thresholding betas based on CI widths.']);
beta = beta(2:end); % remove beta for intercept
beta_thresh = beta; % copy beta vector for thresholding
my_bcis_l = ci(1,2:end); % lower CI
my_bcis_u = ci(2,2:end); % upper CI
sign_l = sign(my_bcis_l); % lower CI signs
sign_u = sign(my_bcis_u); % upper CI signs
sign_dif = sign_l-sign_u; % difference in CI signs
cross_zero = find(sign_dif); % identify zero-crossing CIs
beta_thresh(cross_zero)=NaN; % remove betas with zero-crossing CIs
beta_thresh(beta_thresh == 0)=NaN; % NAN betas with zero values

%%%%%%%%%%%%%%%%%%%%% Outputs %%%%%%%%%%%%%%%%%%%%
% save results to a structure 
disp(['Finished analysis.']);
plsr_results.X_trim = X_trim; % trimmed predictor matrix
plsr_results.X_ind = X_inc_ind; % indices of included cells (i.e. to map back to original predictor matrix)
plsr_results.X_scores = XS; % PLS component scores for retained components
plsr_results.opt_k = opt_k; % optimimum number of PLS components
plsr_results.beta = beta; % beta weights (unthresholded)
plsr_results.beta_thresh = beta_thresh; % beta weights (thresholded)
plsr_results.pred_y = yhat; % predicted values for response variable
plsr_results.r2 = r2; % model R2
plsr_results.r2_ci = ci(:,1); % CIs for R2 value
plsr_results.rng_seed = rng_seed; % rng seed (for reproducibility)
plsr_results.cfg = cfg; % append cfg file for reproducibility

end
