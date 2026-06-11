% This function projects voxel-level weight vectors back to the brain
% and applies thresholding and rescaling if desired.
% Joseph Griffis, 2017, Washington University in St Louis
function [temp] = get_lmwv(temp_nii, lm_has_dmg, lm_weights, rescale_flag)

if isempty(temp_nii)==1
    temp_nii = load_nii('/data/nil-bluearc/corbetta/Studies/SurfaceStroke/Analysis/griffisj/Projects/Structural_FC_Prediction/Lesion_Masks_Fixed_Nov/FCS_026_A_lesion_111_fnirt_333MNI.nii');
end

temp = temp_nii;
temp.img(:,:,:)=NaN;
temp.img(lm_has_dmg) = lm_weights;

if rescale_flag == 1
    temp.img(lm_has_dmg) = temp.img(lm_has_dmg)./max(abs(temp.img(lm_has_dmg)));
    temp.img = temp.img.*10;
elseif rescale_flag == 2
    temp.img(lm_has_dmg) = zscore(temp.img(lm_has_dmg));
elseif rescale_flag == 3
    temp.img(lm_has_dmg) = temp.img(lm_has_dmg).*10000;
else
    
end

end