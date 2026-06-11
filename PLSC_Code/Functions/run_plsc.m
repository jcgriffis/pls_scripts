function plsc = run_plsc(cfg)
% This function computes the PLSC between a predictor matrix and a response
% matrix.

%%%%%%% Inputs:

%%%% cfg structure with the following fields:

% Y: Response matrix
% X: Predictor matrix

%%%%%%% Outputs: 

%%%% plsc structure with the following fields

%%% Fields containing "raw" group-level PLSC results:
% Lx: latent variable scores for X
% Ly: latent variable scores for Y
% Fx: loadings for X
% Fy: loadings for Y
% cve: covariance explained
% U: left singular vectors of R
% S: singular values of R
% V: right singular vectors of R
% nX: rescaled X
% nY: rescaled Y

%%% Fields containing results from permutation testing
% p_vals: permutation p-values for each LV
%%% Fields containing results from bootstrapping
% boot_r: bootstrapped correlations between X and Y scores (rows) for each LV (column)
% bs_ratio_X: bootstrap ratios for each loading (row) on each LV (column)

% Joseph Griffis 2018 Washington University in St Louis

% Variable check
disp('Checking variables...');
if numel(cfg.X(isnan(cfg.X))) > 0
    error('Error: X matrix contains NaN values. Please check data and try again.');
end
if numel(cfg.X(isinf(cfg.X))) > 0
    error('Error: X matrix contains Inf values. Please check data and try again.');
end
if numel(cfg.Y(isnan(cfg.Y))) > 0
    error('Error: Y matrix contains NaN values. Please check data and try again.');
end
if numel(cfg.Y(isinf(cfg.Y))) > 0 
    error('Error: Y matrix contains Inf values. Please check data and try again.');
end
if size(cfg.X,1) ~= size(cfg.Y,1)
    error('Error: X and Y matrices have different numbers of rows. Please check data and try again.');
end
disp('No problems detected with input variables. Running analysis...');

n_obs = size(cfg.X,1);
% Rescale variables
if cfg.zscore == 0
    disp('Z-score flag set to 0: Mean-centering variables...');
    % rescale variables
    plsc.nX = cfg.X - (repmat(mean(cfg.X),[n_obs,1]));
    plsc.nY = cfg.Y - (repmat(mean(cfg.Y),[n_obs,1]));
elseif cfg.zscore == 1
    disp('Z-score flag set to 1: Standardizing variables...');
    % rescale variables
    plsc.nX = cfg.X - (repmat(mean(cfg.X),[n_obs,1]))./(repmat(std(cfg.X),[n_obs,1]));
    plsc.nY = cfg.Y - (repmat(mean(cfg.Y),[n_obs,1]))./(repmat(std(cfg.Y),[n_obs,1]));
end

% compute PLSC
disp('Computing PLSC...');
R = plsc.nY'*plsc.nX; % get covariance of X and Y
if size(plsc.nX, 2) > n_obs || size(plsc.nY, 2) > n_obs
    [plsc.U, plsc.S, plsc.V] = paq(R, n_obs); % compute SVD
else
    [plsc.U, plsc.S, plsc.V] = paq(R); % compute SVDz
end
plsc.Lx = plsc.nX*plsc.V; % get latent variable scores
plsc.Ly = plsc.nY*plsc.U; % get latent variable scores
plsc.cve = (plsc.S.^2)./sum(plsc.S.^2); % compute % covariance explained
plsc.S = diag(plsc.S); % get singular values
plsc.Fx = plsc.V*plsc.S; % get latent variable loadings
plsc.Fy = plsc.U*plsc.S; % get latent variable loadings

% Run permutation tests
plsc = run_perm_plsc(plsc,cfg.n_perm);

% Bootstrap CIs for LV weights on significant LVs
sig_comps = find(plsc.p_vals <= 0.05);
[plsc] = run_boot_plsc(plsc, cfg, sig_comps, n_obs);

disp('Finished running PLSC analysis.');

end
