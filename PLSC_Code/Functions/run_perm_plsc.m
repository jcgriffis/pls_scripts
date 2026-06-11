function [plsc] = run_perm_plsc(plsc, n_perm)

% This function computes runs permutatiosn on the PLSC between a predictor matrix and a response
% matrix.

% Inputs:
% nY: rescaled Response matrix
% nX: rescaled Predictor matrix
% S: singular values from original PLSC
% V: right singular vector from original PLSC
% n_obs: number of observations
% n_perm: number of permutations

% Outputs: 
% sv_resamp: sum of squares of permutation V
% su_resamp: sum of squares of permutation U

% Joseph Griffis 2018, Washington University in St Louis

% Permutation test for significance of LVs 
disp('Running permutation tests...');
n_obs = size(plsc.nX,1);

su_perm = zeros(n_perm, n_obs-1);
for i = 1:n_perm
    
    disp(['Permutation test number ' num2str(i) ' of ' num2str(n_perm)]);
    perm_order = randperm(n_obs)'; % permute order of rows in X
    perm_X = plsc.nX(perm_order,:); % create permuted X matrix 
    R_perm = plsc.nY'*perm_X; % get covariance of permuted X and Y
    
    [U_perm, S_perm, V_perm] = paq(R_perm, size(diag(plsc.S))); % SVD
    if length(S_perm) < length(diag(plsc.S))
        U_perm = [U_perm, zeros(size(U_perm,1), length(diag(plsc.S))-length(S_perm))];
        V_perm = [V_perm, zeros(size(V_perm,1), length(diag(plsc.S))-length(S_perm))];
        S_perm = [S_perm; zeros(length(diag(plsc.S))-length(S_perm),1)];
    end
    
    % Procruste analysis to match LVs 
    [N O P] = paq(plsc.V'*V_perm); % SVD for procrustes analysis
    Q = N*P'; % Get transform 
    S_perm = diag(S_perm); % chagne format for transform
    V_hat_perm = V_perm*S_perm*Q; % apply transform to V
    U_hat_perm = U_perm*S_perm*Q; % apply transform to U
    
    for j = 1:size(U_hat_perm,2) % get effects for permutation tests
       su_perm(i,j) = sqrt(nansum(U_hat_perm(:,j).^2)); % sum of squares of permutation U
    end
end

% Get p-values
plsc.p_vals = [];
for i = 1:size(V_hat_perm,2) % for each LV
    plsc.p_vals(i) = (numel(find(su_perm(:,i) >= plsc.S(i,i)))+1)./(n_perm); % get p-value for LV
end
plsc.p_vals = plsc.p_vals';
plsc.p_vals(plsc.p_vals>1) = 1;
end