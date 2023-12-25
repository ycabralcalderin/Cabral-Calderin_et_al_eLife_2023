opMontage = getOptMontage_0071;
allElec = [reshape([opMontage.LH],numel([opMontage.LH]),1);reshape([opMontage.RH],numel([opMontage.RH]),1)]; %pult all electrodes in a 1 column array

[uni,~,idx] = unique(allElec);
countofE = hist(idx,unique(idx));
figure,imagesc(countofE)
colormap autumn
colorbar
xticks([1:16])
xticklabels(uni')
