function [plsr_results] = run_oos_plsr_fun(cfg)
% This function runs out-of-sample PLSR analyses to predict a response variable (e.g. behavior) from an imaging predictor matrix
% and uses nested leave-one-out cross-validation to determine the optimal number of predictor components and predict held-out cases. 
% Joseph Griffis, 2020

%%%%%%% Inputs %%%%%%%%
% X - an observation-by-predictor (e.g. patient-by-voxel, patient-by-edge), predictor matrix (rows are observations, columns are variables).
% Y - a column vector of responses to be used as the dependent variable.
% trim_flag - a flag for trimming predictor matrix to remove columns with insufficient observations (e.g. voxels that are only lesioned in 1 patient).
% freq_thresh - frequency threshold (i.e. how many patients need non-zero data for a given predictor column to keep it; only used if trim_flag == 1)

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
plsr_results = loocv_plsr_oos(X_trim, Y, cfg);

%%%%%%%%%%%%%%%%%%%%% Outputs %%%%%%%%%%%%%%%%%%%%
% save results to a structure 
disp(['Finished analysis.']);
plsr_results.X_trim = X_trim; % trimmed predictor matrix
plsr_results.X_ind = X_inc_ind; % indices of included cells (i.e. to map back to original predictor matrix)
plsr_results.cfg = cfg; % append cfg file for reproducibility

end