% This script is a general template for running in-sample PLSR analyses 
% It saves thresholded beta maps and plots of % variance explained by each model
% Joseph Griffis 2020

% Inputs should be a predictor matrix X where each row is a subject and eachc column is a variable, and an outcome measure Y corresponding to a column vector (although this script could be modified to accomodate a multi-column Y matrix -- see MATLAB documentation on plsregress for more information).
% X and Y must have the same number of rows.
% X and Y should not contain any NaN or Inf values.
% X should be trimmed to only contain columns with sufficient numbers of non-zero observations.

% This script outputs and saves a plsr_results structure with the following fields:

%%% Fields containing "raw" group-level PLSR results:
% X_trim - trimmed X matrix (i.e. columns removed based on cfg.min_obs and cfg.freq_thresh)
% X_ind - indices of retained columns in original X matrix. (i.e. X_original(:,X_ind) = X_trim
% X_scores - predictor component scores from PLSR model
% opt_k - optimal number of PLS components in final model
% beta - raw (unthresholded) beta weights for original variables in X_trim
% beta_thresh - thresholded beta weights for original variables in X_trim
% pred_y - predicted Y ( pred_y = [ones(length(cfg.Y(cfg.include)), 1) X_trim]*beta
% r2 - r-squared for final model
% r2_ci - CIs for model r-squared
% rng_seed - rng seed information
% cfg -- cfg file containing analysis specifications and data

% clear workspace
clc
clear all

% add paths to relevant functions
plsr_dir = pwd;
addpath(genpath(fullfile(plsr_dir, 'Support_Functions')));
addpath(genpath(fullfile(plsr_dir, 'matlab_nifti')));

% define path to output directory
cfg.out_dir = plsr_dir;

% PLSR settings
cfg.n_boot = 1000; % number of bootstraps
cfg.alpha_thresh = .05; % alpha threshold
cfg.fwe_flag = 1; % apply family-wise error correction
cfg.trim_flag = 1; % trim flag = 1 will remove all columns from the predictor matrix with less than freq_thresh non-zero observations
cfg.min_obs = 1; % minimum value for an observation to be considered valid (must be > 0 if trim flag == 1)
cfg.freq_thresh = 4; % trim flag = 1 will remove all columns from the predictor matrix with less than freq_thresh observations with values > min_obs

% load predictor data
cfg.X_file = fullfile(plsr_dir, 'lesion_data.mat');
load(cfg.X_file); % load file containing predictor data
cfg.X_name = 'lesion'; % predictor matrix name
cfg.X = lesion_data; % lesion matrix (assumes that it's vectorized patient-by-voxel)

% load response data
cfg.Y_file = fullfile(plsr_dir, 'lang.mat');
load(cfg.Y_file); % load file containing response data
cfg.Y_name = 'lang'; % response variable name
cfg.Y = lang; % response variable - expects a column vector

% Define inclusion criteria
cfg.include = find(isnan(cfg.Y)==0); % get indices for patients without missing data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% run PLSR models and save results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run plsr model
[plsr_results] = run_in_sample_plsr_fun(cfg);

% save results
disp('Saving results...');
cd(cfg.out_dir);
save(['plsr_results_' date], 'plsr_results', '-v7.3');

% load nifti file (to fill with beta values)
temp_nii = load_nii([plsr_dir, '/nifti_file.nii']);
temp_nii.hdr.dime.datatype = 16; temp_nii.hdr.dime.bitpix = 32; % make sure data will be stored as continuous values.

% map betas back to volume and save
[beta_img] = get_lmwv(temp_nii, plsr_results.X_ind, plsr_results.beta_thresh, 1); % put back in brain space (note, this will need to be adjusted if you use non-voxel data such as parcels).
save_nii(beta_img, ['lesion_betas_' cfg.Y_name '_freq' num2str(cfg.freq_thresh) '_p' num2str(cfg.alpha_thresh) '.nii']); % save

%% Plot model fit and CIs
fig=bar(plsr_results.r2); % bar plot of R^2
fig.FaceColor='k'; % bar facecolor
box off % remove outline
grid on % turn on grid
set(gca, 'LineWidth', 2); % line width
ylabel('Model R^2'); % y label
xlabel('Model'); % x label
title(cfg.Y_name); % plot title
line([1 1], [plsr_results.r2_ci(:,1)], 'color', [1 0 0], 'LineStyle', '-', 'LineWidth', 2);