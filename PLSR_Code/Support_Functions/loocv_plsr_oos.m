function [results] = loocv_plsr_oos(X,Y,cfg)
% This function identifies the optimal lambda value using leave-one-out
% cross-validation 

% Inputs:
% X: predictor matrix (Nxk)
% Y: response vector (Nx1)
% cfg file
% note - input variables must not include NaN values.
% Outputs: 
% y_pred - fitted Y from final model
% y_real - actual Y
% betas - beta matrix
% pred_err - squared prediction error for final model with optimal lambda
% opt_k - optimal components

% by Joseph Griffis, 2017, Washington University in St Louis
[rows_x, ~] = size(X);
pat_ids = 1:1:rows_x; % create patient indices


% lambda optimization of ridge models
for i = 1:length(pat_ids) % loop through patients
    disp(['PLS Optimization Loop: ' num2str(i)]);
    
    % Get test/train indices for outer loop
    test_outer = i; % current test case
    train_outer = setdiff(pat_ids,test_outer); % training set is all patients but the current patient index
    
    % inner loop to get lambda for this test case
    [opt_k(i)] = get_opt_k(X(train_outer,:),Y(train_outer), cfg);
    
    % split data for outer loop 
    x_train = X(train_outer,:); % training set predictors are Xs for training set indices
    x_test = X(test_outer,:); % test case predictors are Xs for test set index 
    y_train = Y(train_outer); % training Y
    y_test = Y(test_outer); % test Y
    
    % fit ridge regression
    npc = opt_k(i);
    [~,~,~,~,BETA(i,:)] = plsregress(x_train,y_train,npc); % fit PLSR model with j components
    y_pred(i,1) = [1 x_test]*BETA(i,:)'; % get predicted response for test case    
    y_real(i,1) = y_test;
end

% store results for output
numerator = nansum((y_real-y_pred).^2); % get squared error
denominator = nansum((y_real-nanmean(y_pred)).^2);
results.mfit = 1 - (numerator/denominator);
results.y_preds = y_pred;
results.y_real = y_real;
results.betas = BETA;
results.opt_k = opt_k;

end