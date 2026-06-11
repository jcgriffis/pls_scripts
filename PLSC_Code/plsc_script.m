%% PLSC input Script
% This script runs partial least squares correlation (PLSC) analyses (aka BPLS, Symmetric PLSR)
% according to the procedures described in
% Misic et al (2016 -- Cerebral Cortex) and Griffis et al., 2019 (Cell Reports).
% PLSC consists of applying SVD to the cross-covariance matrix computed from two data tables to 
% identify linear combinations of the original variables (latent variables)
% from each table that maximally covary together. The statistical significance of each pair of latent variables
% is determined by permutation testing (i.e. by identifying latent variable pairs whose shared variance is reliably
% greater than that obtained under the permutation null). The statistical significance of the variable loadings on
% the latent variables is determined by computing bootstrapped signal-to-noise ratios (BSRs), which are roughly equivalent
% to z-scores. The correlation between scores (i.e. projected observations) for each pair of latent variables should not 
% be interpreted at face value because it is likely to be inflated, as in canonical correlation analysis.
 
% Inputs should be X and Y matrices where each row is a subject and eachc column is a variable.
% X and Y matrices must have the same number of rows, but may have different numbers of columns.
% X and Y matrices should not contain any NaN or Inf values, as they will cause the SVD to error.
% X and Y matrices should be trimmed to only contain columns with sufficient numbers of non-zero observations.

% This script outputs and saves a plsc structure with the following fields:

%%% Fields containing "raw" group-level PLSC results:
% Lx: latent variable scores for X (patient projections onto latent variables of X)
% Ly: latent variable scores for Y (patient projections onto latent variables of Y)
% Fx: loadings for X (variable loadings on latent variables of X)
% Fy: loadings for Y (variable loadings on latent variables of Y)
% cve: proportion of covariance explained by each pair of latent variables
% U: left singular vectors of covariance matrix
% S: singular values of covariance matrix
% V: right singular vectors of covariance matrix
% nX: rescaled X matrix 
% nY: rescaled Y matrix

%%% Fields containing results from permutation testing
% p_vals: permutation p-values for each pair of latent variables

%%% Fields containing results from bootstrapping
% boot_r: bootstrapped correlations between X and Y scores (rows) for each pair of latent variables (column)
% bs_ratio_X: bootstrap ratios for each variable loading (rows) on each latent variable of X (column)
% bs_ratio_Y: bootstrap ratiso for each variable loading (rows) on each latent variable of Y (column)

%%%%%%%%%%%%%%% Set up and run PLSC analysis %%%%%%%%%%%%%%%%

% Load data and initialize variables
clc
clear all
% Add paths to necessary functions
addpath('/data/nil-bluearc/corbetta/Studies/SurfaceStroke/Analysis/griffisj/Projects/Statistical_Analysis_Code/PLSC_Code/Functions');

% Data directory and file
cfg.data_dir = '/data/nil-bluearc/corbetta/Studies/SurfaceStroke/Analysis/griffisj/Projects/Structural_FC_Prediction/Analyses/Datasets';
cfg.data_file = 'fcsnov_dataset_60vert_5CC_15-Jan-2019';
cfg.out_path = '/data/nil-bluearc/corbetta/Studies/SurfaceStroke/Analysis/griffisj/Projects/Statistical_Analysis_Code/PLSC_Code/';
cfg.out_dir = 'Test_Out';
cd(cfg.out_path);
if ~exist(cfg.out_dir)
    mkdir(cfg.out_dir);
end
load(fullfile(cfg.data_dir, cfg.data_file));

%%%%%% Set analysis parameters
cfg.freq_thresh = 3; % damage threshold for inclusion of column in X matrix (i.e. at least this many observations with values above the binarization threshold -- see next line)
cfg.bin_thresh = 0.01; % binarization threshold for computing column frequencies in predictor matrix (i.e. if equal to 0.01, then only values > 0.01 will be considered when computing the column frequencies)
cfg.n_perm = 1000; % number of permutations for permutation test
cfg.n_boot = 1000; % number of bootstraps for bootstrap SEexc_p = union(exc_p, find(dataset.damage.lesion_side == 1));
cfg.zscore = 0; % Flag to determine whether input matrices should be standardized (cfg.zscore = 1) or mean-centered (cfg.zscore = 0). If standardized, PLSC optimizes based on correlation. If mean-centered, it optimizaes based on covariances.

%%%%%% Set up datasets
include = find(dataset.fc.n_frames.pats >= 180); % inclusion criterion (in this case, # of usable resting frames)
% X and Y matrices
cfg.X = dataset.discon.pats.end.og.s2d(include,:); % declare X matrix
cfg.X_name = 'SDC';
cfg.Y = dataset.fc.mats.pats.z.s2d(include,:); % declare Y matrix 
cfg.Y_name = 'FC';

%%%%%% Clean datasets
% Set NaNs and Infs to 0 in both matrices
cfg.Y(isnan(cfg.Y))=0;
cfg.Y(isinf(cfg.Y))=0;
cfg.X(isnan(cfg.X))=0;
cfg.X(isinf(cfg.X))=0;
% Trim disconnection matrix to remove columns with no or very few observations
[cfg.X, cfg.X_inds] = get_has_dmg(cfg.X, cfg.bin_thresh, cfg.freq_thresh);

%%%%%%% Run PLSC analysis
[plsc] = run_plsc(cfg);

%%
%%%%%%% Save results
disp('Saving results...');
cd(cfg.out_path);
if exist(cfg.out_dir)==0
    mkdir(cfg.out_dir);
end
cd(['./' cfg.out_dir]);
save(['plsc_results_' date], 'plsc', '-v7.3');
save(['plsc_cfg_' date], 'cfg', '-v7.3');