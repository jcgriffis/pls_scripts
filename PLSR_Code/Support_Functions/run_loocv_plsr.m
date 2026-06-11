function [opt_k, p_ress] = run_loocv_plsr(X,Y)
% This function identifies the optimal number of PLSR components using
% leave-one-out cross-validation. 
% Inputs:
% X: predictor matrix (Nxk)
% Y: response vector (Nx1)
% note - input variables must not include NaN values.
% Outputs: 
% opt_k: optimal number of components based on PRESS criterion (Abdi 2010)
% p_ress: predicted residual sum of squares

% by Joseph Griffis, 2018, Washington University in St Louis

[rows_x cols_x] = size(X); % get size of predictor matrix 
dv = (Y); % response variable 
iv = (X); % predictor matrix
pat_ids = 1:1:(rows_x); % create patient indices

dif_press = 1; % initialize dif_press at 1 (criterion for stopping while loop)
j=0; % initialize counter (also corresponds to number of components)

while dif_press > 0 % while the loocv loop with j components has smaller PRESS than the last loop with j-1 components

    j = j+1; % increase counter/n_comp each loop
    
    if j < cols_x
        % run LOOCV loop on patients
        for i = 1:rows_x % loop through patients
            test_case = i; % held out case for testing
            train_set = setdiff(pat_ids,test_case); % training set indices
            x_train = iv(train_set,:); % training set predictors
            x_test = iv(test_case,:); % test case predictors
            y_train = dv(train_set); % training set responses
            y_test = dv(test_case); % test set response
            
            [~,~,~,~,BETA] = plsregress(x_train,y_train,j); % fit PLSR model with j components
            y_pred(i,j) = [1 x_test]*BETA; % get predicted response for test case
            error(i,j) = (y_test-y_pred(i,j)).^2; % get prediction error
        end
        p_ress(j) = nansum(error(:,j)); % get LOOCV predicted residual sum ofsquares for model with j components
        
        if j == 1 % if this is the first loop (i.e. 1 component model), keep going
            dif_press = 1;
        else % else, compute the difference in PRESS between current loop and previous loop
            dif_press = p_ress(j-1)-p_ress(j); % get difference between current PRESS and last PRESS
        end
    else
        dif_press = 0;
    end
end

opt_k = find(p_ress==min(p_ress)); % get optimal component number (is equal to j)

end