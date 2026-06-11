function [X_dmg_only, X_has_dmg] = get_has_dmg(X, bin_thresh, d_thresh)

X_bin = X; X_bin(X_bin >= bin_thresh)=1; X_bin(X_bin<1)=0; % binarize damaged vs. not damaged
X_has_dmg = (nansum(X_bin,1)); X_has_dmg = find(X_has_dmg > d_thresh); % identify damaged connections
X_dmg_only = [X(:,X_has_dmg)]; % remove undamaged connections from X

end
