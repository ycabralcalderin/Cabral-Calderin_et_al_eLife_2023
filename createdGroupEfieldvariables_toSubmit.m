MainFolder = uigetdir([], 'choose main directory with all subjects'); %Main directory with all subjects


subjects = {''};

BOLDfolder = 'GroupData2';%folder with fMRI contrast data

results_folderS = fullfile('simulations_FC6_TP8-P8_FC5_TP7-P7', 'mni_volumes'); %folder with standard simulations
results_folderI = fullfile('simulations_Func_0071_V2', 'mni_volumes'); %folder with IndivMontage simulations
results_folderR = fullfile('simulations_T7-T8_ring2', 'mni_volumes'); %folder with IndivMontage simulations


field_name = 'normE';

mni_image_suffix = ['_TDCS_1_scalar_MNI_' field_name '.nii.gz'];

% template_image = nifti_load(fullfile(MainFolder, subjects{1}, results_folder, [subjects{1} mni_image_suffix]));
% field_avg = zeros(size(template_image.vol));
% field_avgalls = zeros([length(subjects) size(template_image.vol)]);

for i = 1:length(subjects)
    
    subj = subjects{i};
    % load the nifti images
    EfieldGroup.imgS(i,:,:,:) = nifti_load(fullfile(MainFolder, subj, results_folderS, [subj mni_image_suffix]));
    EfieldGroup.imgI(i,:,:,:) = nifti_load(fullfile(MainFolder, subj, results_folderI, [subj mni_image_suffix]));
    EfieldGroup.imgR(i,:,:,:) = nifti_load(fullfile(MainFolder, subj, results_folderR, [subj mni_image_suffix]));

end
save EfieldGroup_0071 EfieldGroup -v7.3



