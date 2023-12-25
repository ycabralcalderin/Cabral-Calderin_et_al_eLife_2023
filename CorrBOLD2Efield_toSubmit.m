MainFolder = uigetdir([], 'choose main directory with all subjects'); %Main directory with all subjects


SUBJlist = {'SUBJ1',2,4,4,3; 'SUBJ2',1,4,[],3;'SUBJ3',2,5,2,3;'SUBJ4',2,5,5,3;'SUBJ5',2,3,1,3;'SUBJ6',2,1,4,3;'SUBJ7',1,4,[],3;'SUBJ8',2,5,3,3;
    'SUBJ9',2,1,3,3;'SUBJ10',2,2,1,3;'SUBJ11',2,2,4,3;'SUBJ12',2,4,3,3;'SUBJ13',2,2,2,3;'SUBJ13',2,1,1,1;'SUBJ15',2,5,3,3;'SUBJ16',2,2,3,1;'SUBJ17',2,2,3,1;
    'SUBJ18',2,1,5,1;'SUBJ19',2,1,3,2;'SUBJ20',2,4,3,2;'SUBJ21',2,4,5,2;'SUBJ22',2,5,4,1;'SUBJ23',2,5,1,1;'SUBJ24',2,4,2,1;'SUBJ25',2,3,2,1;'SUBJ26',2,1,2,3;
    'SUBJ27',2,5,1,1;'SUBJ28',2,4,4,1;'SUBJ29',2,1,3,1;'SUBJ30',2,1,3,1;'SUBJ31',2,2,1,2;'SUBJ32',2,2,2,1;'SUBJ33',2,5,5,1;'SUBJ34',2,3,5,1;'SUBJ35',2,2,5,1;
    'SUBJ36',2,1,4,1;'SUBJ37',2,3,3,2;'SUBJ38',2,1,4,2;'SUBJ39',1,1,[],1;'SUBJ40',2,3,5,2;'SUBJ41',2,2,2,1;'SUBJ42',2,4,2,2};

indm       = [1 2 3 4 5 6 7 8 9 10 11 12 13 21 26 38 40 42];%subjects with individualized montage
inds       = [14 15 16 17 18 19 20 22 23 24 25 27 28 29 30 31 32 33 34 35 36 37 39 41];%subjects with standard montage
missingMRI = [15 28 33];%
BOLDfolder = '';%folder with fMRI contrast data

results_folderS = fullfile('simulations_FC6_TP8-P8_FC5_TP7-P7', 'mni_volumes'); %folder with standard simulations
results_folderI = fullfile('simulations_Func_0071_V2', 'mni_volumes'); %folder with IndivMontage simulations
results_folderR = fullfile('simulations_T7-T8_ring2', 'mni_volumes'); %folder with IndivMontage simulations
load('GroupEfieldAna_orig_wholeHem_0071_V2.mat')


field_name = 'normE';

mni_image_suffix = ['_TDCS_1_scalar_MNI_' field_name '.nii.gz'];

% template_image = nifti_load(fullfile(MainFolder, subjects{1}, results_folder, [subjects{1} mni_image_suffix]));
% field_avg = zeros(size(template_image.vol));
% field_avgalls = zeros([length(subjects) size(template_image.vol)]);

for i = 1:length(subjects)
    
    subj = subjects{i};
    % load the nifti images
    imgS = nifti_load(fullfile(MainFolder, subj, results_folderS, [subj mni_image_suffix]));
    imgI = nifti_load(fullfile(MainFolder, subj, results_folderI, [subj mni_image_suffix]));
    imgR = nifti_load(fullfile(MainFolder, subj, results_folderR, [subj mni_image_suffix]));

    boldDIR = dir(fullfile(BOLDfolder, [subj '_con_*rs.nii.gz']));
    imgBOLD = nifti_load(fullfile(boldDIR.folder,boldDIR.name));
    funcMaskDIR = dir(fullfile(BOLDfolder, [subj '_AudiovsBase_rs.nii.gz']));
    imgFunMask = nifti_load(fullfile(funcMaskDIR.folder,funcMaskDIR.name));
    maskGM = nifti_load(fullfile(MainFolder, subj,['m2m_' subj '/toMNI/gm_MNI.nii.gz']));
    maskWM = nifti_load(fullfile(MainFolder, subj,['m2m_' subj '/toMNI/wm_MNI.nii.gz']));
  
%create a new mask by finding the voxels where values are not nan for
%anatomical masks, BOLD contrast and Efield
nanWM = maskWM.vol>0;
nanGM = maskGM.vol>0;
nanAnat = nanGM+nanWM; 
nanAnat(nanAnat>0)=1;%mask include volxels were either WM or GM are true
nanBOLD = ~isnan(imgBOLD.vol); %find missing values in contrast(nan)
nanS = ~isnan(imgS.vol);%find missing values in standard montage(nan)
nanI = ~isnan(imgI.vol);%find missing values in individualized montage(nan)
comMask = nanAnat+nanBOLD+nanS+nanI; %combine all masks
comMask(comMask<4)=0;
comMask(comMask==4)=1; %final mask will be 1 if all images have nonnan values, otherwise it will be cero
maskBOLDpos = imgBOLD.vol>0;
maskBOLDneg = imgBOLD.vol<=0;
%some plotting
figure
subplot(2,3,1)
imagesc(squeeze(imgS.vol(:,111,:))')
title('standard montage')
subplot(2,3,2)
imagesc(squeeze(imgI.vol(:,111,:))')
title('individual montage')
subplot(2,3,3)
imagesc(squeeze(imgBOLD.vol(:,111,:))')
title('BOLD')
subplot(2,3,4)
imagesc(squeeze(maskGM.vol(:,111,:))')
title('GM')
subplot(2,3,5)
imagesc(squeeze(maskWM.vol(:,111,:))')
title('WM')
subplot(2,3,6)
imagesc(squeeze(imgFunMask.vol(:,111,:))')
title('FuncMask')
figure
subplot(2,2,1)
imagesc(squeeze(comMask(:,111,:))')
title('final Mask')
subplot(2,2,2)
imagesc(squeeze(nanAnat(:,111,:))')
title('Anat Mask')
subplot(2,2,3)
imagesc(squeeze(nanBOLD(:,111,:))')
title('BOLD Mask')
subplot(2,2,4)
imagesc(squeeze(nanS(:,111,:))')
title('S Mask')

comMask = logical(comMask);

GroupEfieldall(i).I = imgI.vol;
GroupEfieldall(i).S = imgS.vol;
GroupEfieldall(i).R = imgR.vol;
GroupEfieldall(i).B = imgBOLD.vol;

GroupCorr(i,1) = corr(imgI.vol(comMask),imgS.vol(comMask)); %correlation individual to standard
GroupCorr(i,2) = corr(imgI.vol(comMask),imgBOLD.vol(comMask));%individual to bold
GroupCorr(i,3) = corr(imgS.vol(comMask),imgBOLD.vol(comMask));%standard to bold
GroupCorr(i,4) = corr(imgI.vol(comMask),imgR.vol(comMask));
GroupCorr(i,5) = corr(imgS.vol(comMask),imgR.vol(comMask));
GroupCorr(i,6) = corr(imgR.vol(comMask),imgBOLD.vol(comMask));

corrMatrix(i,:,:) = [corr(imgI.vol(comMask),imgI.vol(comMask)),corr(imgI.vol(comMask),imgS.vol(comMask)),corr(imgI.vol(comMask),imgR.vol(comMask));...
    corr(imgS.vol(comMask),imgI.vol(comMask)),corr(imgS.vol(comMask),imgS.vol(comMask)),corr(imgS.vol(comMask),imgR.vol(comMask));...
    corr(imgR.vol(comMask),imgI.vol(comMask)),corr(imgR.vol(comMask),imgS.vol(comMask)),corr(imgR.vol(comMask),imgR.vol(comMask))];
funcMask = imgFunMask.vol;
funcMask(funcMask>0)=1;
%extract data from functional mask
GroupEfieldFuncROI(i,1) = mean(imgI.vol(logical(funcMask))); 
GroupEfieldFuncROI(i,2) = mean(imgS.vol(logical(funcMask))); 
GroupEfieldFuncROI(i,3) = mean(imgR.vol(logical(funcMask))); 

GroupEfieldFuncpos(i,1) = mean(imgI.vol(logical(maskBOLDpos))); 
GroupEfieldFuncpos(i,2) = mean(imgS.vol(logical(maskBOLDpos))); 
GroupEfieldFuncpos(i,3) = mean(imgR.vol(logical(maskBOLDpos))); 

GroupEfieldFuncneg(i,1) = mean(imgI.vol(logical(maskBOLDneg))); 
GroupEfieldFuncneg(i,2) = mean(imgS.vol(logical(maskBOLDneg))); 
GroupEfieldFuncneg(i,3) = mean(imgR.vol(logical(maskBOLDneg))); 

%end
end

%correlatetion across participants
for s1 =1:length(subjects)
    GroupEfieldI(s1,:,:,:)=GroupEfieldall(s1).I;
    GroupEfieldS(s1,:,:,:)=GroupEfieldall(s1).S;
    GroupEfieldR(s1,:,:,:)=GroupEfieldall(s1).R;
    GroupEfieldB(s1,:,:,:)=GroupEfieldall(s1).B;

    for s2 =1:length(subjects)
        corrMI(s1,s2)=corr(reshape(GroupEfieldall(s1).I,numel(GroupEfieldall(s1).I),1),reshape(GroupEfieldall(s2).I,numel(GroupEfieldall(s2).I),1), 'rows','complete');
        corrMS(s1,s2)=corr(reshape(GroupEfieldall(s1).S,numel(GroupEfieldall(s1).S),1),reshape(GroupEfieldall(s2).S,numel(GroupEfieldall(s2).S),1), 'rows','complete');
        corrMR(s1,s2)=corr(reshape(GroupEfieldall(s1).R,numel(GroupEfieldall(s1).R),1),reshape(GroupEfieldall(s2).R,numel(GroupEfieldall(s2).R),1), 'rows','complete');
        corrMB(s1,s2)=corr(reshape(GroupEfieldall(s1).B,numel(GroupEfieldall(s1).B),1),reshape(GroupEfieldall(s2).B,numel(GroupEfieldall(s2).B),1), 'rows','complete');
    end
end


figure
subplot(2,2,1)
imagesc(corrMI)
caxis([-1 1])
colorbar
title('spatial correlation Indiv Montage')
subplot(2,2,2)
imagesc(corrMS)
caxis([-1 1])

colorbar
title('spatial correlation Standard Montage')
subplot(2,2,3)
imagesc(corrMR)
caxis([-1 1])

colorbar
title('spatial correlation HD Montage')
subplot(2,2,4)
imagesc(corrMB)
caxis([-1 1])

colorbar
title('spatial correlation BOLD')
mean(corrMI,'all')
mean(corrMS,'all')
mean(corrMR,'all')
mean(corrMB,'all')

figure, subplot(1,2,1)
imagesc(GroupCorr)
xlabel('Condition')
xticklabels({'I2S','I2B','S2B','I2R','S2R','R2B'})
ylabel('Subjects')
colorbar
title('correlation Efield2BOLD')
subplot(1,2,2)
imagesc(GroupCorr(:,2:6))
xlabel('Condition')
ylabel('Subjects')
colorbar
title('correlation Efield2BOLD')
xticklabels({'I2B','S2B','I2R','S2R','R2B'})

EfieldBOLDdata.corr = GroupCorr;
EfieldBOLDdata.corrMatrix = corrMatrix;
EfieldBOLDdata.Read = {'I2S','I2BOLD','S2BOLD' 'I2R' 'S2R' 'R2B'};
GroupCorr(i,1) = corr(imgI.vol(comMask),imgS.vol(comMask)); %correlation individual to standard
GroupCorr(i,2) = corr(imgI.vol(comMask),imgBOLD.vol(comMask));%individual to bold
GroupCorr(i,3) = corr(imgS.vol(comMask),imgBOLD.vol(comMask));%standard to bold
GroupCorr(i,4) = corr(imgI.vol(comMask),imgR.vol(comMask));
GroupCorr(i,5) = corr(imgS.vol(comMask),imgR.vol(comMask));
GroupCorr(i,6) = corr(imgR.vol(comMask),imgBOLD.vol(comMask));

save(fullfile(MainFolder,'EFIELD2BOLD_correlation'),'EfieldBOLDdata')
save(fullfile(MainFolder,'EFIELD2BOLD_GroupEfieldFuncROI'),'GroupEfieldFuncROI')
save(fullfile(MainFolder,'EFIELD2BOLD_GroupEfieldFuncpos'),'GroupEfieldFuncpos')
save(fullfile(MainFolder,'EFIELD2BOLD_GroupEfieldFuncneg'),'GroupEfieldFuncneg')
GroupCorrZ = atanh(EfieldBOLDdata.corr);
corrMatrixZ = atanh(EfieldBOLDdata.corrMatrix);

figure 
distributionPlotYCC([GroupCorrZ(:,2:3) GroupCorrZ(:,6)],'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
%xticklabels({'I2S';'I2B';'S2B';'I2R';'S2R';'R2B'})
yticklabels({'I2B';'S2B';'R2B'})
xlim([-0.3 0.5])
title('Correlation Efield-to-BOLD activation')
ylabel('Montage')
xlabel('Fishers z-scores')



[statsAll.CORR_anovaP,statsAll.CORR_anovaTable]= anova_rm([GroupCorrZ(:,2:3) GroupCorrZ(:,6)]);
[statsAll.CorrI2S_H2,statsAll.CorrI2S_P2,statsAll.CorrI2S_CI2,statsAll.CorrI2S_STATS2] = ttest(GroupCorrZ(:,2),GroupCorrZ(:,3));
statsAll.CorrI2S_pBonf =statsAll.CorrI2S_P2*3;
[statsAll.CorrI2R_H2,statsAll.CorrI2R_P2,statsAll.CorrI2R_CI2,statsAll.CorrI2R_STATS2] = ttest(GroupCorrZ(:,2),GroupCorrZ(:,6));
statsAll.CorrI2R_pBonf =statsAll.CorrI2R_P2*3;
[statsAll.CorrS2R_H2,statsAll.CorrS2R_P2,statsAll.CorrS2R_CI2,statsAll.CorrS2R_STATS2] = ttest(GroupCorrZ(:,3),GroupCorrZ(:,6));
statsAll.CorrS2R_pBonf =statsAll.CorrS2R_P2*3;

squeeze(mean(corrMatrix))
mean2plot=squeeze(mean(corrMatrixZ));
mean2plot(mean2plot>3)=nan;%change infinitiy for nan
figure
imagesc(mean2plot)
caxis([0 1.5])
colormap(copper)
colorbar
xticklabels({'I';'S';'R'})
yticklabels({'I';'S';'R'})

cmap = colormap(copper);
% [H2,P2,CI2,STATS2] = ttest(GroupEfieldFuncROI(:,1),GroupEfieldFuncROI(:,2));
% [H3,P3,CI3,STATS3] = ttest(GroupEfieldFuncpos(:,1),GroupEfieldFuncpos(:,2));
% [H4,P4,CI4,STATS4] = ttest(GroupEfieldFuncneg(:,1),GroupEfieldFuncneg(:,2));
[statsAll.ROI_anovaP,statsAll.ROI_anovaTable] = anova_rm(GroupEfieldFuncROI);
[statsAll.I2S_H2,statsAll.I2S_P2,statsAll.I2S_CI2,statsAll.I2S_STATS2] = ttest(GroupEfieldFuncROI(:,1),GroupEfieldFuncROI(:,2));
statsAll.I2S_pBonf =statsAll.I2S_P2*3;
[statsAll.S2R_H2,statsAll.S2R_P2,statsAll.S2R_CI2,statsAll.S2R_STATS2] = ttest(GroupEfieldFuncROI(:,2),GroupEfieldFuncROI(:,3));
statsAll.S2R_pBonf =statsAll.S2R_P2*3;
[statsAll.I2R_H2,statsAll.I2R_P2,statsAll.I2R_CI2,statsAll.I2R_STATS2] = ttest(GroupEfieldFuncROI(:,1),GroupEfieldFuncROI(:,3));
statsAll.I2R_pBonf =statsAll.I2R_P2*3;

figure, 
subplot(1,3,1)
distributionPlotYCC(GroupEfieldFuncROI,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
hold on
%plot(GroupEfieldFuncROI')
yticklabels({'Individual';'Standard';'Ring'})
title('funcROI')
xlim([-0.05 0.35])
subplot(1,3,2)
distributionPlotYCC(GroupEfieldFuncpos,'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
hold on
%plot(GroupEfieldFuncROI')
yticklabels({'Individual';'Standard';'Ring'})
title('funcPOS')
subplot(1,3,3)
distributionPlotYCC(GroupEfieldFuncneg,'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
hold on
%plot(GroupEfieldFuncROI')
yticklabels({'Individual';'Standard';'Ring'})
title('funcNeg')
save(fullfile(MainFolder,'EFIELD2BOLD_Stats'),'statsAll')

% for x=1:39
% figure
% subplot(2,2,1)
% imagesc(squeeze(GroupEfieldall(x).I(:,100,:)))
% subplot(2,2,2),imagesc(squeeze(GroupEfieldall(x).S(:,100,:)))
% subplot(2,2,3),imagesc(squeeze(GroupEfieldall(x).R(:,100,:)))
% subplot(2,2,4),imagesc(squeeze(GroupEfieldall(x).B(:,100,:)))
% end
%1.- jacknife get distance from indiv to group mean
for s1 =1:length(subjects)
subjallid     = ones(length(subjects),1);
subjallid(s1) = 0;
BOLDmean               = mean(GroupEfieldB(logical(subjallid),:,:,:),'omitnan');
BOLDSubjdist(s1,:,:,:) = GroupEfieldB(s1,:,:,:)-BOLDmean;
[BOLDSubjcorr2Mean(s1,1),BOLDSubjcorr2Mean(s1,2)] = corr(reshape(GroupEfieldB(s1,:,:,:),1,numel(GroupEfieldB(s1,:,:,:)))',reshape(BOLDmean,1,numel(BOLDmean))', 'rows','complete');
Imean                  = mean(GroupEfieldI(logical(subjallid),:,:,:),'omitnan');
ISubjdist(s1,:,:,:)    = GroupEfieldI(s1,:,:,:)-Imean;
[ISubjcorr2Mean(s1,1),ISubjcorr2Mean(s1,2)] = corr(reshape(GroupEfieldI(s1,:,:,:),1,numel(GroupEfieldI(s1,:,:,:)))',reshape(Imean,1,numel(Imean))', 'rows','complete');
Smean                  = mean(GroupEfieldS(logical(subjallid),:,:,:),'omitnan');
SSubjdist(s1,:,:,:)    = GroupEfieldS(s1,:,:,:)-Smean;
[SSubjcorr2Mean(s1,1),SSubjcorr2Mean(s1,2)] = corr(reshape(GroupEfieldS(s1,:,:,:),1,numel(GroupEfieldS(s1,:,:,:)))',reshape(Smean,1,numel(Smean))', 'rows','complete');
Rmean                  = mean(GroupEfieldR(logical(subjallid),:,:,:),'omitnan');
RSubjdist(s1,:,:,:)    = GroupEfieldR(s1,:,:,:)-Rmean;
[RSubjcorr2Mean(s1,1),RSubjcorr2Mean(s1,2)] = corr(reshape(GroupEfieldR(s1,:,:,:),1,numel(GroupEfieldR(s1,:,:,:)))',reshape(Rmean,1,numel(Rmean))', 'rows','complete');
end
Sub2GroupCorrZ = atanh([BOLDSubjcorr2Mean(:,1) ISubjcorr2Mean(:,1) SSubjcorr2Mean(:,1) RSubjcorr2Mean(:,1)]);
for s1 =1:length(subjects)
BOLDSubjmean(s1) = mean(abs(BOLDSubjdist(s1,:,:,:)),'all','omitnan');
ISubjmean(s1) = mean(abs(ISubjdist(s1,:,:,:)),'all','omitnan');
SSubjmean(s1) = mean(abs(SSubjdist(s1,:,:,:)),'all','omitnan');
RSubjmean(s1) = mean(abs(RSubjdist(s1,:,:,:)),'all','omitnan');
end
[r,p]=corr(BOLDSubjmean',ISubjmean','type','Spearman');
figure
subplot(2,3,1)
scatter(BOLDSubjmean,ISubjmean)
subplot(2,3,2)
scatter(BOLDSubjmean,SSubjmean)
subplot(2,3,3)
scatter(BOLDSubjmean,RSubjmean)
subplot(2,3,4)
scatter(ISubjmean,SSubjmean)
subplot(2,3,5)
scatter(ISubjmean,RSubjmean)
subplot(2,3,6)
scatter(SSubjmean,RSubjmean)
figure
subplot(2,3,1)
scatter(Sub2GroupCorrZ(:,1),Sub2GroupCorrZ(:,2))
subplot(2,3,2)
scatter(Sub2GroupCorrZ(:,1),Sub2GroupCorrZ(:,3))
subplot(2,3,3)
scatter(Sub2GroupCorrZ(:,1),Sub2GroupCorrZ(:,4))
subplot(2,3,4)
scatter(Sub2GroupCorrZ(:,2),Sub2GroupCorrZ(:,3))
subplot(2,3,5)
scatter(Sub2GroupCorrZ(:,2),Sub2GroupCorrZ(:,4))
subplot(2,3,6)
scatter(Sub2GroupCorrZ(:,3),Sub2GroupCorrZ(:,4))

%%%%%%USING THE SPECIFIC MONATGE USED IN THE TACS session
indm = [1 2 3 4 5 6 7 8 9 10 11 12 13 20 25 35 37 39];%subjects with individualized montage

inds = [14 15 16 17 18 19 21 22 23 24 26 27 28 29 30 31 32 33 34 36 38];%subjects with standard montage

%2.---compute distance between Efield peak and BOLD ROI center coordinate
load('ROIcentroid')
for i = indm %1:length(ROIcentroid) individual montage group
if ~strcmp(subjects{i},ROIcentroid(i).ID)
error 'WARNING...SUBJECT ID DOES NOT MATCH'
else
if ROIcentroid(i).roi(1,1)<0 %left ROI
templfcoord = ROIcentroid(i).roi(1,:);
temprfcoord = ROIcentroid(i).roi(2,:);
else
templfcoord = ROIcentroid(i).roi(2,:);
temprfcoord = ROIcentroid(i).roi(1,:);
end
Efield2BOLDdist(i,1) = sqrt((templfcoord(1) -FieldVal.EnormLH(i,1).peak_pos(1,1))^2 + (templfcoord(2) -FieldVal.EnormLH(i,1).peak_pos(1,2))^2+(templfcoord(3) -FieldVal.EnormLH(i,1).peak_pos(1,3))^2); %compute euclidian distance
Efield2BOLDdist(i,2) = sqrt((temprfcoord(1) -FieldVal.EnormRH(i,1).peak_pos(1,1))^2 + (temprfcoord(2) -FieldVal.EnormRH(i,1).peak_pos(1,2))^2+(temprfcoord(3) -FieldVal.EnormRH(i,1).peak_pos(1,3))^2); %compute euclidian distance
Efield2BOLDdist(i,3) = mean(Efield2BOLDdist(i,1:2));
Efield2BOLDdist(i,4) = 1;
end
end
for i = inds %1:length(ROIcentroid) STANDARD montage group
if ~strcmp(subjects{i},ROIcentroid(i).ID)
error 'WARNING...SUBJECT ID DOES NOT MATCH'
else
if ROIcentroid(i).roi(1,1)<0 %left ROI
templfcoord = ROIcentroid(i).roi(1,:);
temprfcoord = ROIcentroid(i).roi(2,:);
else
templfcoord = ROIcentroid(i).roi(2,:);
temprfcoord = ROIcentroid(i).roi(1,:);
end
Efield2BOLDdist(i,1) = sqrt((templfcoord(1) -FieldVal.EnormLH(i,2).peak_pos(1,1))^2 + (templfcoord(2) -FieldVal.EnormLH(i,2).peak_pos(1,2))^2+(templfcoord(3) -FieldVal.EnormLH(i,2).peak_pos(1,3))^2); %compute euclidian distance
Efield2BOLDdist(i,2) = sqrt((temprfcoord(1) -FieldVal.EnormRH(i,2).peak_pos(1,1))^2 + (temprfcoord(2) -FieldVal.EnormRH(i,2).peak_pos(1,2))^2+(temprfcoord(3) -FieldVal.EnormRH(i,2).peak_pos(1,3))^2); %compute euclidian distance
Efield2BOLDdist(i,3) = mean(Efield2BOLDdist(i,1:2));
Efield2BOLDdist(i,4) = 2;
end
end

%distance to target for each montage
for i = 1:length(ROIcentroid)
if ~strcmp(subjects{i},ROIcentroid(i).ID)
error 'WARNING...SUBJECT ID DOES NOT MATCH'
else
if ROIcentroid(i).roi(1,1)<0 %left ROI
templfcoord = ROIcentroid(i).roi(1,:);
temprfcoord = ROIcentroid(i).roi(2,:);
else
templfcoord = ROIcentroid(i).roi(2,:);
temprfcoord = ROIcentroid(i).roi(1,:);
end
Efield2BOLDdistI(i,1) = sqrt((templfcoord(1) -FieldVal.EnormLH(i,1).peak_pos(1,1))^2 + (templfcoord(2) -FieldVal.EnormLH(i,1).peak_pos(1,2))^2+(templfcoord(3) -FieldVal.EnormLH(i,1).peak_pos(1,3))^2); %compute euclidian distance
Efield2BOLDdistI(i,2) = sqrt((temprfcoord(1) -FieldVal.EnormRH(i,1).peak_pos(1,1))^2 + (temprfcoord(2) -FieldVal.EnormRH(i,1).peak_pos(1,2))^2+(temprfcoord(3) -FieldVal.EnormRH(i,1).peak_pos(1,3))^2); %compute euclidian distance
Efield2BOLDdistI(i,3) = mean(Efield2BOLDdistI(i,1:2));
Efield2BOLDdistI(i,4) = 1;

Efield2BOLDdistS(i,1) = sqrt((templfcoord(1) -FieldVal.EnormLH(i,2).peak_pos(1,1))^2 + (templfcoord(2) -FieldVal.EnormLH(i,2).peak_pos(1,2))^2+(templfcoord(3) -FieldVal.EnormLH(i,2).peak_pos(1,3))^2); %compute euclidian distance
Efield2BOLDdistS(i,2) = sqrt((temprfcoord(1) -FieldVal.EnormRH(i,2).peak_pos(1,1))^2 + (temprfcoord(2) -FieldVal.EnormRH(i,2).peak_pos(1,2))^2+(temprfcoord(3) -FieldVal.EnormRH(i,2).peak_pos(1,3))^2); %compute euclidian distance
Efield2BOLDdistS(i,3) = mean(Efield2BOLDdistS(i,1:2));
Efield2BOLDdistS(i,4) = 1;

Efield2BOLDdistR(i,1) = sqrt((templfcoord(1) -FieldVal.EnormLH(i,3).peak_pos(1,1))^2 + (templfcoord(2) -FieldVal.EnormLH(i,3).peak_pos(1,2))^2+(templfcoord(3) -FieldVal.EnormLH(i,3).peak_pos(1,3))^2); %compute euclidian distance
Efield2BOLDdistR(i,2) = sqrt((temprfcoord(1) -FieldVal.EnormRH(i,3).peak_pos(1,1))^2 + (temprfcoord(2) -FieldVal.EnormRH(i,3).peak_pos(1,2))^2+(temprfcoord(3) -FieldVal.EnormRH(i,3).peak_pos(1,3))^2); %compute euclidian distance
Efield2BOLDdistR(i,3) = mean(Efield2BOLDdistR(i,1:2));
Efield2BOLDdistR(i,4) = 1;
end
end

figure,
subplot(1,2,1)
distributionPlotYCC([Efield2BOLDdistI(:,3) Efield2BOLDdistS(:,3) Efield2BOLDdistR(:,3)],'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
title('Distance PeakEfield2Target')
xlim([0 25])
ylabel({'I';'S'; 'R'})
subplot(2,2,2)
distributionPlotYCC(Efield2BOLDdist(Efield2BOLDdist(:,4)==1,3),'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
title('Distance PeakEfield2Target Individual Group')
xlim([0 25])
subplot(2,2,4)
distributionPlotYCC(Efield2BOLDdist(Efield2BOLDdist(:,4)==2,3),'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
title('Distance PeakEfield2Target Standard Group')
xlim([0 25])

[statsAll.EfieldBoldDistAnovap,statsAll.EfieldBoldDistAnovatable]=anova_rm([Efield2BOLDdistI(:,3),Efield2BOLDdistS(:,3),Efield2BOLDdistR(:,3)]);
[statsAll.EfieldBoldDisth(1),statsAll.EfieldBoldDistp(1),statsAll.EfieldBoldDistci(1,:),statsAll.EfieldBoldDiststats(1).stats]=ttest(Efield2BOLDdistI(:,3),Efield2BOLDdistS(:,3));
[statsAll.EfieldBoldDisth(2),statsAll.EfieldBoldDistp(2),statsAll.EfieldBoldDistci(2,:),statsAll.EfieldBoldDiststats(2).stats]=ttest(Efield2BOLDdistI(:,3),Efield2BOLDdistR(:,3));
[statsAll.EfieldBoldDisth(3),statsAll.EfieldBoldDistp(3),statsAll.EfieldBoldDistci(3,:),statsAll.EfieldBoldDiststats(3).stats]=ttest(Efield2BOLDdistS(:,3),Efield2BOLDdistR(:,3));
statsAll.EfieldBoldDistpBonf(1)=statsAll.EfieldBoldDistp(1)*3;
statsAll.EfieldBoldDistpBonf(2)=statsAll.EfieldBoldDistp(2)*3;
statsAll.EfieldBoldDistpBonf(3)=statsAll.EfieldBoldDistp(3)*3;

%3.---%COMPUTINg OPTIMALITY INDEX
%square root of the sum of squares of intensity, focality, BOLD overlap,
%distance BOL peak 2 Efield peak
%convert values to z-score before any computation
for s = 1:length(subjects)
EnormG(s,:)=[FieldVal.Enorm(s,1).p95 FieldVal.Enorm(s,2).p95];
EnormAL(s,:)=[FieldVal.Enormal(s,1).p95 FieldVal.Enormal(s,2).p95];
EnormF(s,:)=[FieldVal.FOC(s,1).c50 FieldVal.FOC(s,2).c50];
Corr2BOLD(s,:)=EfieldBOLDdata.corr(s,2:3);
end
EnormGnormalize = (EnormG - min(min(EnormG)))/(max(max(EnormG)) - min(min(EnormG)));
EnormalGnormalize = (EnormAL - min(min(EnormAL)))/(max(max(EnormAL)) - min(min(EnormAL)));

%first invert FOC by substracting from the max
EnormF2 = -EnormF;%multiply for -1 to make it positivity correlated with better
EnormFnormalize = (EnormF2 - min(min(EnormF2)))/(max(max(EnormF2)) - min(min(EnormF2)));
EnormROInormalize = (GroupEfieldFuncROI - min(min(GroupEfieldFuncROI)))/(max(max(GroupEfieldFuncROI)) - min(min(GroupEfieldFuncROI)));
%normalize(GroupEfieldFuncROI(:,1:2),'range');
%Corr2BOLDfisherZ = atanh(Corr2BOLD);
Corr2BOLDrange1 =(Corr2BOLD-min(min(Corr2BOLD))./(max(max(Corr2BOLD))-min(min(Corr2BOLD))));
Efield2BOLDdist2 = -Efield2BOLDdist;%multiply for -1 to make it positivity correlated with better
Efield2BOLDdistnormalize = (Efield2BOLDdist2 - min(min(Efield2BOLDdist2)))/(max(max(Efield2BOLDdist2)) - min(min(Efield2BOLDdist2)));
%normalize(Efield2BOLDdist2(:,3));
% z=.5.*log((1+Corr2BOLD)./(1-Corr2BOLD)); %other way to coput fisher
%
% Corr2BOLDfisherZ-z
for i = indm %1:length(ROIcentroid) STANDARD montage group
Optimality1(i,1) = Efield2BOLDdistnormalize(i) + EnormROInormalize(i,1) + EnormFnormalize(i,1); %compute euclidian distance
Optimality2(i,1) = Efield2BOLDdistnormalize(i) + EnormGnormalize(i,1) + EnormFnormalize(i,1); %compute euclidian distance
Optimality3(i,1) = Corr2BOLDrange1(i,1) + EnormROInormalize(i,1) + EnormFnormalize(i,1); %compute euclidian distance
Optimality4(i,1) = Corr2BOLDrange1(i,1) + EnormGnormalize(i,1) + EnormFnormalize(i,1); %compute euclidian distance
Optimality5(i,1) = Efield2BOLDdistnormalize(i) + EnormalGnormalize(i,1) + EnormFnormalize(i,1); %compute euclidian distance
Optimality6(i,1) = Corr2BOLDrange1(i,1) + EnormalGnormalize(i,1) + EnormFnormalize(i,1); %compute euclidian distance
OtherVariables(i,:) = [Corr2BOLDrange1(i,1) EnormROInormalize(i,1) EnormGnormalize(i,1) EnormFnormalize(i,1) EnormalGnormalize(i,1)];
end
for i = inds %1:length(ROIcentroid) STANDARD montage group
Optimality1(i,1) = Efield2BOLDdistnormalize(i) + EnormROInormalize(i,2) + EnormFnormalize(i,2); %compute euclidian distance
Optimality2(i,1) = Efield2BOLDdistnormalize(i) + EnormGnormalize(i,2) + EnormFnormalize(i,2); %compute euclidian distance
Optimality3(i,1) = Corr2BOLDrange1(i,2) + EnormROInormalize(i,2) + EnormFnormalize(i,2); %compute euclidian distance
Optimality4(i,1) = Corr2BOLDrange1(i,2) + EnormGnormalize(i,2) + EnormFnormalize(i,2); %compute euclidian distance
Optimality5(i,1) = Efield2BOLDdistnormalize(i) + EnormalGnormalize(i,2) + EnormFnormalize(i,2); %compute euclidian distance
Optimality6(i,1) = Corr2BOLDrange1(i,2) + EnormalGnormalize(i,2) + EnormFnormalize(i,2); %compute euclidian distance
OtherVariables(i,:) = [Corr2BOLDrange1(i,2) EnormROInormalize(i,2) EnormGnormalize(i,2) EnormFnormalize(i,2) EnormalGnormalize(i,2)];
end
figure,
subplot(2,6,1)
title('Optimality1: Dist-EnormROI-Foc')
distributionPlotYCC(Optimality1(Efield2BOLDdist(:,4)==1),'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])
subplot(2,6,7)
distributionPlotYCC(Optimality1(Efield2BOLDdist(:,4)==2),'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])
subplot(2,6,2)
title('Optimality2: Dist-Enorm-Foc')
distributionPlotYCC(Optimality2(Efield2BOLDdist(:,4)==1),'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])
subplot(2,6,8)
distributionPlotYCC(Optimality2(Efield2BOLDdist(:,4)==2),'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])
subplot(2,6,3)
title('Optimality3: CorrBold-EnormROI-Foc')
distributionPlotYCC(Optimality3(Efield2BOLDdist(:,4)==1),'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])
subplot(2,6,9)
distributionPlotYCC(Optimality3(Efield2BOLDdist(:,4)==2),'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])
subplot(2,6,4)
title('Optimality4: CorrBold-Enorm-Foc')
distributionPlotYCC(Optimality4(Efield2BOLDdist(:,4)==1),'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])
subplot(2,6,10)
distributionPlotYCC(Optimality4(Efield2BOLDdist(:,4)==2),'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])
subplot(2,6,5)
title('Optimality5: Dist-Enormal-Foc')
distributionPlotYCC(Optimality5(Efield2BOLDdist(:,4)==1),'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])
subplot(2,6,11)
distributionPlotYCC(Optimality5(Efield2BOLDdist(:,4)==2),'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])
subplot(2,6,6)
title('Optimality6: CorrBold-Enormal-Foc')
distributionPlotYCC(Optimality6(Efield2BOLDdist(:,4)==1),'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])
subplot(2,6,12)
distributionPlotYCC(Optimality6(Efield2BOLDdist(:,4)==2),'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])

[statsAll.hop1,statsAll.pop1,statsAll.ciop1,statsAll.statsop1]=ttest2(Optimality1(Efield2BOLDdist(:,4)==1),Optimality1(Efield2BOLDdist(:,4)==2));
[statsAll.hop2,statsAll.pop2,statsAll.ciop2,statsAll.statsop2]=ttest2(Optimality2(Efield2BOLDdist(:,4)==1),Optimality2(Efield2BOLDdist(:,4)==2));
[statsAll.hop3,statsAll.pop3,statsAll.ciop3,statsAll.statsop3]=ttest2(Optimality3(Efield2BOLDdist(:,4)==1),Optimality3(Efield2BOLDdist(:,4)==2));
[statsAll.hop4,statsAll.pop4,statsAll.ciop4,statsAll.statsop4]=ttest2(Optimality4(Efield2BOLDdist(:,4)==1),Optimality4(Efield2BOLDdist(:,4)==2));
[statsAll.hop5,statsAll.pop5,statsAll.ciop5,statsAll.statsop5]=ttest2(Optimality5(Efield2BOLDdist(:,4)==1),Optimality5(Efield2BOLDdist(:,4)==2));
[statsAll.hop6,statsAll.pop6,statsAll.ciop6,statsAll.statsop6]=ttest2(Optimality6(Efield2BOLDdist(:,4)==1),Optimality6(Efield2BOLDdist(:,4)==2));

%%%%optimality for all montages
for s = 1:length(subjects)
EnormGA(s,:)=[FieldVal.Enorm(s,1).p95 FieldVal.Enorm(s,2).p95 FieldVal.Enorm(s,3).p95];
EnormALGA(s,:)=[FieldVal.Enormal(s,1).p95 FieldVal.Enormal(s,2).p95 FieldVal.Enormal(s,3).p95];
EnormFA(s,:)=[FieldVal.FOC(s,1).c50 FieldVal.FOC(s,2).c50 FieldVal.FOC(s,3).c50];
Corr2BOLDA(s,:)=EfieldBOLDdata.corr(s,[2:3 6]);
end
EnormGnormalizeA = (EnormGA - min(min(EnormGA)))/(max(max(EnormGA)) - min(min(EnormGA)));%normalize(EnormGA,'range');%normalize from 0-1 for the whole matrix
EnormalGnormalizeA = (EnormALGA - min(min(EnormALGA)))/(max(max(EnormALGA)) - min(min(EnormALGA)));

%first invert FOC by substracting from the max
EnormF2A = -EnormFA;%multiply for -1 to make it positivity correlated with better
EnormFnormalizeA = (EnormF2A - min(min(EnormF2A)))/(max(max(EnormF2A)) - min(min(EnormF2A)));
EnormROInormalizeA = (GroupEfieldFuncROI - min(min(GroupEfieldFuncROI)))/(max(max(GroupEfieldFuncROI)) - min(min(GroupEfieldFuncROI)));%normalize(GroupEfieldFuncROI,'range');
%Corr2BOLDfisherZA = atanh(Corr2BOLDA);
Corr2BOLDrangeA =(Corr2BOLDA-min(min(Corr2BOLDA))./(max(max(Corr2BOLDA))-min(min(Corr2BOLDA))));
Efield2BOLDdistt =[Efield2BOLDdistI(:,3) Efield2BOLDdistS(:,3) Efield2BOLDdistR(:,3)];
Efield2BOLDdist2A = -Efield2BOLDdistt;%multiply for -1 to make it positivity correlated with better
Efield2BOLDdistnormalizeA = (Efield2BOLDdist2A - min(min(Efield2BOLDdist2A)))./(max(max(Efield2BOLDdist2A)) - min(min(Efield2BOLDdist2A)));%normalize(GroupEfieldFuncROI,'range');

figure 
subplot(2,3,1),distributionPlotYCC(Efield2BOLDdistnormalizeA,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
title('Efield to BOLD distance')
xlabel('distance(a.u.)')
ylabel('Montage')
yticklabels({'I', 'S', 'R'})
subplot(2,3,2),distributionPlotYCC(Corr2BOLDrangeA,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
title('Efield to BOLD correlation')
xlabel('correlation(fishers zscores)')
ylabel('Montage')
yticklabels({'I', 'S', 'R'})
subplot(2,3,3),distributionPlotYCC(EnormROInormalizeA,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
title('Enorm ROI')
xlabel('|E|(a.u.)')
ylabel('Montage')
yticklabels({'I', 'S', 'R'})
subplot(2,3,4),distributionPlotYCC(EnormGnormalizeA,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
title('Enorm all')
xlabel('|E|(a.u.)')
ylabel('Montage')
yticklabels({'I', 'S', 'R'})
subplot(2,3,5),distributionPlotYCC(EnormFnormalizeA,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
title('Focality')
xlabel('Focality(a.u.)')
ylabel('Montage')
yticklabels({'I', 'S', 'R'})

figure 
subplot(2,3,1),distributionPlotYCC(Efield2BOLDdistt,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
title('Efield to BOLD distance')
xlabel('distance(mm)')
ylabel('Montage')
yticklabels({'I', 'S', 'R'})
xlim([0 25])
subplot(2,3,2),distributionPlotYCC(Corr2BOLDA,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
title('Efield to BOLD correlation')
xlabel('correlation(fishers zscores)')
ylabel('Montage')
yticklabels({'I', 'S', 'R'})
xlim([-0.3 0.45])
subplot(2,3,3),distributionPlotYCC(GroupEfieldFuncROI,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
title('Enorm ROI')
xlabel('|E|(V/m)')
ylabel('Montage')
yticklabels({'I', 'S', 'R'})
xlim([0 0.4])

subplot(2,3,4),distributionPlotYCC(EnormGA,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
title('Enorm all')
xlabel('|E|(V/m)')
ylabel('Montage')
yticklabels({'I', 'S', 'R'})
xlim([0 0.2])

subplot(2,3,5),distributionPlotYCC(EnormALGA,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
title('Enormal all')
xlabel('E(V/m)')
ylabel('Montage')
yticklabels({'I', 'S', 'R'})
xlim([0 0.1])

subplot(2,3,6),distributionPlotYCC(EnormFA,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
title('Focality')
xlabel('Focality(mm)')
ylabel('Montage')
yticklabels({'I', 'S', 'R'})
xlim([2000 130000])

[statsAll.EnormGA_anovaP,statsAll.EnormGA_anovaTable]= anova_rm(EnormGA);
[statsAll.EnormGAI2S_H2,statsAll.EnormGAI2S_P2,statsAll.EnormGAI2S_CI2,statsAll.EnormGAI2S_STATS2] = ttest(EnormGA(:,1),EnormGA(:,2));
statsAll.EnormGAI2S_pBonf =statsAll.EnormGAI2S_P2*3;
[statsAll.EnormGAI2R_H2,statsAll.EnormGAI2R_P2,statsAll.EnormGAI2R_CI2,statsAll.EnormGAI2R_STATS2] = ttest(EnormGA(:,1),EnormGA(:,3));
statsAll.EnormGAI2R_pBonf =statsAll.EnormGAI2R_P2*3;
[statsAll.EnormGAS2R_H2,statsAll.EnormGAS2R_P2,statsAll.EnormGAS2R_CI2,statsAll.EnormGAS2R_STATS2] = ttest(EnormGA(:,2),EnormGA(:,3));
statsAll.EnormGAS2R_pBonf =statsAll.EnormGAS2R_P2*3;

[statsAll.EnormALGA_anovaP,statsAll.EnormALGA_anovaTable]= anova_rm(EnormALGA);
[statsAll.EnormALGAI2S_H2,statsAll.EnormALGAI2S_P2,statsAll.EnormALGAI2S_CI2,statsAll.EnormALGAI2S_STATS2] = ttest(EnormALGA(:,1),EnormALGA(:,2));
statsAll.EnormALGAI2S_pBonf =statsAll.EnormALGAI2S_P2*3;
[statsAll.EnormALGAI2R_H2,statsAll.EnormALGAI2R_P2,statsAll.EnormALGAI2R_CI2,statsAll.EnormALGAI2R_STATS2] = ttest(EnormALGA(:,1),EnormALGA(:,3));
statsAll.EnormALGAI2R_pBonf =statsAll.EnormALGAI2R_P2*3;
[statsAll.EnormALGAS2R_H2,statsAll.EnormALGAS2R_P2,statsAll.EnormALGAS2R_CI2,statsAll.EnormALGAS2R_STATS2] = ttest(EnormALGA(:,2),EnormALGA(:,3));
statsAll.EnormALGAS2R_pBonf =statsAll.EnormALGAS2R_P2*3;

[statsAll.FOC_anovaP,statsAll.FOC_anovaTable]= anova_rm(EnormFA);
[statsAll.FOCI2S_H2,statsAll.FOCI2S_P2,statsAll.FOCI2S_CI2,statsAll.FOCI2S_STATS2] = ttest(EnormFA(:,1),EnormFA(:,2));
statsAll.FOCI2S_pBonf =statsAll.FOCI2S_P2*3;
[statsAll.FOCI2R_H2,statsAll.FOCI2R_P2,statsAll.FOCI2R_CI2,statsAll.FOCI2R_STATS2] = ttest(EnormFA(:,1),EnormFA(:,3));
statsAll.FOCI2R_pBonf =statsAll.FOCI2R_P2*3;
[statsAll.FOCS2R_H2,statsAll.FOCS2R_P2,statsAll.FOCS2R_CI2,statsAll.FOCS2R_STATS2] = ttest(EnormFA(:,2),EnormFA(:,3));
statsAll.FOCS2R_pBonf =statsAll.FOCS2R_P2*3;

%FIT repeated measures model
within = table([1;2;3],'VariableNames',{'montage'});
between = table(Optimality1,Optimality2,Optimality3,...
'VariableNames',{'IO','S','HD'});
rm = fitrm(between,'IO-HD ~ 1','WithinDesign',within);
ranovatbl = ranova(rm);
anova_rm([Optimality1,Optimality2,Optimality3])

for i = 1:length(ROIcentroid) 
    for m=1:3
Optimality1A(i,m) = Efield2BOLDdistnormalizeA(i,m) + EnormROInormalizeA(i,m) + EnormFnormalizeA(i,m); %compute euclidian distance
Optimality2A(i,m) = Efield2BOLDdistnormalizeA(i,m) + EnormGnormalizeA(i,m) + EnormFnormalizeA(i,m); %compute euclidian distance
Optimality3A(i,m) = Corr2BOLDrangeA(i,m) + EnormROInormalizeA(i,m) + EnormFnormalizeA(i,m); %compute euclidian distance
Optimality4A(i,m) = Corr2BOLDrangeA(i,m) + EnormGnormalizeA(i,m) + EnormFnormalizeA(i,m); %compute euclidian distance
Optimality5A(i,m) = Efield2BOLDdistnormalizeA(i,m) + EnormalGnormalizeA(i,m) + EnormFnormalizeA(i,m); %compute euclidian distance
Optimality6A(i,m) = Corr2BOLDrangeA(i,m) + EnormalGnormalizeA(i,m) + EnormFnormalizeA(i,m); %compute euclidian distance
    end
end
figure, 
subplot(2,3,1),distributionPlotYCC(Optimality1A,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])
title('Optimality1: Dist-EnormROI-Foc')
ylabel('Montage')
yticklabels({'I', 'S', 'R'})
subplot(2,3,2),distributionPlotYCC(Optimality2A,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])
title('Optimality2: Dist-Enorm-Foc')
ylabel('Montage')
yticklabels({'I', 'S', 'R'})
subplot(2,3,3),distributionPlotYCC(Optimality3A,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])
title('Optimality3: Corr-EnormROI-Foc')
ylabel('Montage')
yticklabels({'I', 'S', 'R'})
subplot(2,3,4),distributionPlotYCC(Optimality4A,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])
title('Optimality4: Corr-Enorm-Foc')
ylabel('Montage')
yticklabels({'I', 'S', 'R'})
subplot(2,3,5),distributionPlotYCC(Optimality5A,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])
title('Optimality5: Dist-Enormal-Foc')
ylabel('Montage')
yticklabels({'I', 'S', 'R'})
subplot(2,3,6),distributionPlotYCC(Optimality6A,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])
title('Optimality6: Corr-Enormal-Foc')
ylabel('Montage')
yticklabels({'I', 'S', 'R'})


[statsAll.Optimality1Anovap,statsAll.Optimality1Anovatable]=anova_rm(Optimality1A);
[statsAll.Optimality1h(1),statsAll.Optimality1p(1),statsAll.Optimality1ci(1,:),statsAll.Optimality1stats(1).stats]=ttest(Optimality1A(:,1),Optimality1A(:,2));
[statsAll.Optimality1h(2),statsAll.Optimality1p(2),statsAll.Optimality1ci(2,:),statsAll.Optimality1stats(2).stats]=ttest(Optimality1A(:,1),Optimality1A(:,3));
[statsAll.Optimality1h(3),statsAll.Optimality1p(3),statsAll.Optimality1ci(3,:),statsAll.Optimality1stats(3).stats]=ttest(Optimality1A(:,2),Optimality1A(:,3));
statsAll.Optimality1pBonf(1)=statsAll.Optimality1p(1)*3;
statsAll.Optimality1pBonf(2)=statsAll.Optimality1p(2)*3;
statsAll.Optimality1pBonf(3)=statsAll.Optimality1p(3)*3;

[statsAll.Optimality2Anovap,statsAll.Optimality2Anovatable]=anova_rm(Optimality2A);
[statsAll.Optimality2h(1),statsAll.Optimality2p(1),statsAll.Optimality2ci(1,:),statsAll.Optimality2stats(1).stats]=ttest(Optimality2A(:,1),Optimality2A(:,2));
[statsAll.Optimality2h(2),statsAll.Optimality2p(2),statsAll.Optimality2ci(2,:),statsAll.Optimality2stats(2).stats]=ttest(Optimality2A(:,1),Optimality2A(:,3));
[statsAll.Optimality2h(3),statsAll.Optimality2p(3),statsAll.Optimality2ci(3,:),statsAll.Optimality2stats(3).stats]=ttest(Optimality2A(:,2),Optimality2A(:,3));
statsAll.Optimality2pBonf(1)=statsAll.Optimality2p(1)*3;
statsAll.Optimality2pBonf(2)=statsAll.Optimality2p(2)*3;
statsAll.Optimality2pBonf(3)=statsAll.Optimality2p(3)*3;

[statsAll.Optimality3Anovap,statsAll.Optimality3Anovatable]=anova_rm(Optimality3A);
[statsAll.Optimality3h(1),statsAll.Optimality3p(1),statsAll.Optimality3ci(1,:),statsAll.Optimality3stats(1).stats]=ttest(Optimality3A(:,1),Optimality3A(:,2));
[statsAll.Optimality3h(2),statsAll.Optimality3p(2),statsAll.Optimality3ci(2,:),statsAll.Optimality3stats(2).stats]=ttest(Optimality3A(:,1),Optimality3A(:,3));
[statsAll.Optimality3h(3),statsAll.Optimality3p(3),statsAll.Optimality3ci(3,:),statsAll.Optimality3stats(3).stats]=ttest(Optimality3A(:,2),Optimality3A(:,3));
statsAll.Optimality3pBonf(1)=statsAll.Optimality3p(1)*3;
statsAll.Optimality3pBonf(2)=statsAll.Optimality3p(2)*3;
statsAll.Optimality3pBonf(3)=statsAll.Optimality3p(3)*3;

[statsAll.Optimality4Anovap,statsAll.Optimality4Anovatable]=anova_rm(Optimality4A);
[statsAll.Optimality4h(1),statsAll.Optimality4p(1),statsAll.Optimality4ci(1,:),statsAll.Optimality4stats(1).stats]=ttest(Optimality4A(:,1),Optimality4A(:,2));
[statsAll.Optimality4h(2),statsAll.Optimality4p(2),statsAll.Optimality4ci(2,:),statsAll.Optimality4stats(2).stats]=ttest(Optimality4A(:,1),Optimality4A(:,3));
[statsAll.Optimality4h(3),statsAll.Optimality4p(3),statsAll.Optimality4ci(3,:),statsAll.Optimality4stats(3).stats]=ttest(Optimality4A(:,2),Optimality4A(:,3));
statsAll.Optimality4pBonf(1)=statsAll.Optimality4p(1)*3;
statsAll.Optimality4pBonf(2)=statsAll.Optimality4p(2)*3;
statsAll.Optimality4pBonf(3)=statsAll.Optimality4p(3)*3;

[statsAll.Optimality5Anovap,statsAll.Optimality5Anovatable]=anova_rm(Optimality5A);
[statsAll.Optimality5h(1),statsAll.Optimality5p(1),statsAll.Optimality5ci(1,:),statsAll.Optimality5stats(1).stats]=ttest(Optimality5A(:,1),Optimality5A(:,2));
[statsAll.Optimality5h(2),statsAll.Optimality5p(2),statsAll.Optimality5ci(2,:),statsAll.Optimality5stats(2).stats]=ttest(Optimality5A(:,1),Optimality5A(:,3));
[statsAll.Optimality5h(3),statsAll.Optimality5p(3),statsAll.Optimality5ci(3,:),statsAll.Optimality5stats(3).stats]=ttest(Optimality5A(:,2),Optimality5A(:,3));
statsAll.Optimality5pBonf(1)=statsAll.Optimality5p(1)*3;
statsAll.Optimality5pBonf(2)=statsAll.Optimality5p(2)*3;
statsAll.Optimality5pBonf(3)=statsAll.Optimality5p(3)*3;

[statsAll.Optimality6Anovap,statsAll.Optimality6Anovatable]=anova_rm(Optimality6A);
[statsAll.Optimality6h(1),statsAll.Optimality6p(1),statsAll.Optimality6ci(1,:),statsAll.Optimality6stats(1).stats]=ttest(Optimality6A(:,1),Optimality6A(:,2));
[statsAll.Optimality6h(2),statsAll.Optimality6p(2),statsAll.Optimality6ci(2,:),statsAll.Optimality6stats(2).stats]=ttest(Optimality6A(:,1),Optimality6A(:,3));
[statsAll.Optimality6h(3),statsAll.Optimality6p(3),statsAll.Optimality6ci(3,:),statsAll.Optimality6stats(3).stats]=ttest(Optimality6A(:,2),Optimality6A(:,3));
statsAll.Optimality6pBonf(1)=statsAll.Optimality6p(1)*3;
statsAll.Optimality6pBonf(2)=statsAll.Optimality6p(2)*3;
statsAll.Optimality6pBonf(3)=statsAll.Optimality6p(3)*3;

[r,p]=corr(BOLDSubjmean',Optimality1,'type','Spearman');
[r,p]=corr(Optimality1,SSubjmean','type','Spearman');
[r,p]=corr(BOLDSubjmean',Optimality1A(:,1)-Optimality1A(:,2),'type','Spearman');
[r,p]=corr(atanh(BOLDSubjcorr2Mean(:,1)),Optimality1A(:,1)-Optimality1A(:,2),'type','Spearman');

figure,scatter(ones(size(Optimality1(Efield2BOLDdist(:,4)==1))),Optimality1(Efield2BOLDdist(:,4)==1))
hold on
scatter(ones(size(Optimality1(Efield2BOLDdist(:,4)==2))),Optimality1(Efield2BOLDdist(:,4)==2))


% figure,
% subplot(2,3,1)
% distributionPlotYCC(Optimality1(:,1),'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
% hold on
% distributionPlotYCC(Optimality1(:,2),'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
% distributionPlotYCC(Optimality1(:,3),'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
% title('Optimality1')
% xlim([0 4])
% subplot(2,3,2)
% distributionPlotYCC(Optimality2(:,1),'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
% hold on
% distributionPlotYCC(Optimality2(:,2),'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
% distributionPlotYCC(Optimality2(:,3),'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
% title('Optimality1')
% xlim([0 4])
% subplot(2,3,3)
% distributionPlotYCC(Optimality3(:,1),'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
% hold on
% distributionPlotYCC(Optimality3(:,2),'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
% distributionPlotYCC(Optimality3(:,3),'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
% title('Optimality1')
% xlim([0 4])
% subplot(2,3,4)
% distributionPlotYCC(Optimality4(:,1),'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
% hold on
% distributionPlotYCC(Optimality4(:,2),'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
% distributionPlotYCC(Optimality4(:,3),'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
% title('Optimality1')
% xlim([0 4])
% subplot(2,3,5)
% distributionPlotYCC(Optimality5(:,1),'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
% hold on
% distributionPlotYCC(Optimality5(:,2),'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
% distributionPlotYCC(Optimality5(:,3),'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
% title('Optimality1')
% xlim([0 4])
% subplot(2,3,6)
% distributionPlotYCC(Optimality6(:,1),'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
% hold on
% distributionPlotYCC(Optimality6(:,2),'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
% distributionPlotYCC(Optimality6(:,3),'showMM',0,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
% title('Optimality1')
% xlim([0 4])
for subj =1:length(FieldVal.EnormR)
    for m =1:size(FieldVal.EnormR,2)
    normERatio(subj,m) = FieldVal.EnormR(subj,m).p95./FieldVal.EnormL(subj,m).p95;
    normalERatio(subj,m) = FieldVal.EnormalR(subj,m).p95./FieldVal.EnormalL(subj,m).p95;

    end
end
figure, distributionPlotYCC(normERatio,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])
title('normE ratio')
ylabel('Montage')
yticklabels({'I', 'S', 'R'})

figure, distributionPlotYCC(normalERatio,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
xlim([0.5 3])
title('normal E ratio')
ylabel('Montage')
yticklabels({'I', 'S', 'R'})

%create new variable with important values for tACS analysis
save('statsAll')

save('EfieldVar2tACS_0071_V2','Optimality1','Optimality2','Optimality3','Optimality4','Optimality5','Optimality6','OtherVariables','Efield2BOLDdist','Sub2GroupCorrZ')