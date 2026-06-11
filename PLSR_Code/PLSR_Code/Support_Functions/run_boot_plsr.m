% This function serves as input to the bootci anonymous function.
% Outputs are a matrix of CIs for model fit and betas. 
% Beta CIs don't include the
% intercept.

% Joseph Griffis, 2018, Washington University in St. Louis    

function [out] = run_boot_plsr(X,Y,k)
    [~,~,~,~,BETA] = plsregress(X,Y,k); % fit fixed effects model with opt_k   
    yhat = [ones(length(Y), 1) X]*BETA; % get fitted Y
    numerator = nansum((Y-yhat).^2); % get squared error
    denominator = nansum((Y-nanmean(Y)).^2);
    gof = 1 - (numerator/denominator);
    betas = BETA(2:end); % get non-intercept betas
    out = [gof;  betas]; % output matrix
end