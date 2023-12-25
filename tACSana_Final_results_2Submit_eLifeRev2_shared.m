%pipeline for analysing FM-stimulus+tACS data
clear
close all
clc

%DEFINE DIRECTORIES
MAINdir  = uigetdir (); %main path to raw data, e.g., '/Users/yuranny.cabral/Documents/gapDetectionProject/tACS/DATA'
BEHdir   = fullfile(MAINdir, 'Beh'); %raw beh data
EEGdir   = fullfile(MAINdir, 'tACS'); %raw eeg data (for tACS signal)
OUTdir   = fullfile(MAINdir, 'ANA'); %main output directory
cd(OUTdir)
% 
load('SUBJlist')

indm = [1 2 3 4 5 6 7 8 9 10 11 12 13 21 26 38 40 42];%subjects with individualized montage
inds = [14 15 16 17 18 19 20 22 23 24 25 27 28 29 30 31 32 33 34 35 36 37 39 41];%subjects with standard montage
missingMRI = [15 28 33];

%define some general variables
NrGapBins    = 9; %number of FM stimulus phase bins for gap presentation
NrAudtACSlag = 6; %tACS-Audio phase lag

ft_defaults
NrIte = 1000;
%load data with the info about detection FM phase and tACS phase instead of
%raw values
load (fullfile(OUTdir,'dataGroupallGaps.mat'),'dataGroupallGaps')

%preallocate some variables
dataGroupPool = struct([]);
data = struct([]);

%run the analysis using the chosen #of bins and percent overlap
TACSlagwidth = 0; %no overlap between bins

%parameters for binning and for fitting the cosine
phasebinGap   = 0:2*pi/NrGapBins:2*pi;
lagBins       = 0:2*pi/NrAudtACSlag:2*pi; %0:2*pi/4:2*pi; %this is the phase in the center of the bin
maxSeparation = 2*pi/NrAudtACSlag/2; %this is  the half width of a bin for not overlapping
lagBins1      =  wrapTo2Pi(lagBins-maxSeparation-TACSlagwidth/2);
lagBins2      =  wrapTo2Pi(lagBins+maxSeparation+TACSlagwidth/2);

X         = 0:2*pi/NrGapBins:2*pi-2*pi/NrGapBins; %this is my phase vector for fitting the cosine functions
X         = X+2*pi/NrGapBins/2;
lag       = pi;
intercept = 0.5;
amp       = 0.5;
params    = [lag intercept amp];
          
for subj=1:length(dataGroupallGaps)%loop across subjects
       close all
       if isempty(dataGroupallGaps(2).firstShamTr)
           NrSess         = 1;
       else
           NrSess         = 2;
       end
    gapDetect      = [];
    gaptACSangle   = [];
    gapAUDangle    = [];
    gapAUD_tACSlag = [];
    gapStimcond    = [];
    
    for Session = 1:NrSess
        if ~isempty(dataGroupallGaps(subj,Session).gapDetect)
           
            gapDetect      = dataGroupallGaps(subj,Session).gapDetect;
            gaptACSangle   = dataGroupallGaps(subj,Session).gaptACSangle';
            gapAUDangle    = dataGroupallGaps(subj,Session).gapAUDangle;
            gapAUD_tACSlag = dataGroupallGaps(subj,Session).gapAUD_tACSlag';
            gapStimcond    = dataGroupallGaps(subj,Session).gapStimcond;
            
            %first separate sham from real
            %main
            tACS1angleGapMain  = wrapTo2Pi(gaptACSangle(gapStimcond==1)); %tACS phase at gap onset
            trCurrGapTrackMain = gapDetect(gapStimcond==1); %gap info from eeg triggers
            gapAUDtACSlagMain  = wrapTo2Pi(gapAUD_tACSlag(gapStimcond==1)); %phase lag between tACS and audio signal
            gapAUDphaseMain    = gapAUDangle(gapStimcond==1); %FM phase at gap onset
            %sham
            tACS1angleGapsham  = wrapTo2Pi(gaptACSangle(gapStimcond==2)); %tACS phase at gap onset
            trCurrGapTracksham = gapDetect(gapStimcond==2); %gap info from eeg triggers
            gapAUDtACSlagsham  = wrapTo2Pi(gapAUD_tACSlag(gapStimcond==2)); %phase lag between tACS and audio signal
            gapAUDphasesham    = gapAUDangle(gapStimcond==2); %FM phase at gap onset
            
            for gapBin = 1:length(phasebinGap)-1
                currGapbintACS  = find(tACS1angleGapMain>=phasebinGap(gapBin)&tACS1angleGapMain<phasebinGap(gapBin+1))'; %gaps in current tACS phase bin
                currGapbinStim  = find(gapAUDphaseMain>=phasebinGap(gapBin)&gapAUDphaseMain<phasebinGap(gapBin+1))'; %gaps im current AUD bin
                               length(currGapbintACS)
                               length(currGapbinStim)
                for b=1:length(lagBins)-1 %separate by tACS-AUD phase lag
                    if lagBins1(b)>lagBins2(b)
                        currGaptACSbinSlag = find((tACS1angleGapMain>=phasebinGap(gapBin)&tACS1angleGapMain<phasebinGap(gapBin+1))&(gapAUDtACSlagMain>=lagBins1(b)|gapAUDtACSlagMain<lagBins2(b)))';
                        currGapStimbinSlag = find((gapAUDphaseMain>=phasebinGap(gapBin)&gapAUDphaseMain<phasebinGap(gapBin+1))&(gapAUDtACSlagMain>=lagBins1(b)|gapAUDtACSlagMain<lagBins2(b)))';
                    else
                        currGaptACSbinSlag = find((tACS1angleGapMain>=phasebinGap(gapBin)&tACS1angleGapMain<phasebinGap(gapBin+1))&...
                            (gapAUDtACSlagMain>=lagBins1(b)&gapAUDtACSlagMain<lagBins2(b)))';
                        currGapStimbinSlag =find((gapAUDphaseMain>=phasebinGap(gapBin)&gapAUDphaseMain<phasebinGap(gapBin+1))&...
                            (gapAUDtACSlagMain>=lagBins1(b)& gapAUDtACSlagMain<lagBins2(b)))';
                    end
                                        length(currGaptACSbinSlag)
                                        length(currGapStimbinSlag)
                    dataGroupPool(subj,Session).hrbytACSbinStimLagMain(b,gapBin) = mean(trCurrGapTrackMain(currGaptACSbinSlag));
                    dataGroupPool(subj,Session).hrbyStimbintACSLagMain(b,gapBin) = mean(trCurrGapTrackMain(currGapStimbinSlag));
                    
                                        figure, subplot(2,2,1), circ_plot(tACS1angleGapMain(currGaptACSbinSlag))
                                        title(['TACS gap bin: ' num2str(gapBin) ' lag bin' num2str(b)])
                                        subplot(2,2,2), circ_plot(gapAUDtACSlagMain(currGaptACSbinSlag))
                                        title('lag')
                                        subplot(2,2,3), circ_plot(gapAUDphaseMain(currGapStimbinSlag))
                                        title(['AUDIO gap bin: ' num2str(gapBin) ' lag bin' num2str(b)])
                                        subplot(2,2,4), circ_plot(gapAUDtACSlagMain(currGapStimbinSlag))
                                        title('lag')
                end
                %separate by both FM and tACS signal
                for gapBin2 = 1:length(phasebinGap)-1
                    currGapbinAT = find((tACS1angleGapMain>=phasebinGap(gapBin)&tACS1angleGapMain<phasebinGap(gapBin+1))&...
                        (gapAUDphaseMain>=phasebinGap(gapBin2)&gapAUDphaseMain<phasebinGap(gapBin2+1)))';
                    
                    %  length(currGapbinAT)
                    dataGroupPool(subj,Session).hrbyStimbytACSbinMain(gapBin2,gapBin) = mean(trCurrGapTrackMain(currGapbinAT));
                    dataGroupPool(subj,Session).hrbyStimbytACSbinMainNr(gapBin2,gapBin) = length(currGapbinAT);
                end
                dataGroupPool(subj,Session).hrbyStimbinMain(gapBin) = mean(trCurrGapTrackMain(currGapbinStim));
                dataGroupPool(subj,Session).hrbytACSbinMain(gapBin) = mean(trCurrGapTrackMain(currGapbintACS));
            end
            
            for gapBin = 1:length(phasebinGap)-1
                
                currGapbintACSsham = find(tACS1angleGapsham>=phasebinGap(gapBin)&tACS1angleGapsham<phasebinGap(gapBin+1))';
                currGapbinStimsham = find(gapAUDphasesham>=phasebinGap(gapBin)&gapAUDphasesham<phasebinGap(gapBin+1))';
                                length(currGapbintACSsham)
                                length(currGapbinStimsham)
                dataGroupPool(subj,Session).hrbyStimbinsham(gapBin) = mean(trCurrGapTracksham(currGapbinStimsham));
                dataGroupPool(subj,Session).hrbytACSbinsham(gapBin) = mean(trCurrGapTracksham(currGapbintACSsham));
                
                %separate by both FM and tACS signal
                for gapBin2 = 1:length(phasebinGap)-1
                    currGapbinATs = find((tACS1angleGapsham>=phasebinGap(gapBin)&tACS1angleGapsham<phasebinGap(gapBin+1))&...
                        (gapAUDphasesham>=phasebinGap(gapBin2)&gapAUDphasesham<phasebinGap(gapBin2+1)))';
                    %  length(currGapbinATs)
                    dataGroupPool(subj,Session).hrbyStimbytACSbinsham(gapBin2,gapBin) = mean(trCurrGapTracksham(currGapbinATs));
                    dataGroupPool(subj,Session).hrbyStimbytACSbinshamNr(gapBin2,gapBin) = length(currGapbinATs);
                    
                end
            end
            for l = 1:length(lagBins)-1
                try
                    %fit by gap phase and stim_tacs phase lag
                    [dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT(l,1:3), dataGroupPool(subj,Session).hrbytACSbinStimLagMain_RESNORM(l,1),dataGroupPool(subj,Session).hrbytACSbinStimLagMain_RESIDUAL(l,1:18),dataGroupPool(subj,Session).hrbytACSbinStimLagMain_prefPhase(l,1)] = fitCos2RatebyPhase ([dataGroupPool(subj,Session).hrbytACSbinStimLagMain(l,:) dataGroupPool(subj,Session).hrbytACSbinStimLagMain(l,:)]',X,params);
                    dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT_yAll(l,1:9) = dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT(l,2) + dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT(l,3).*(cos(X + dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT(l,1))); % if you want to plot the predicted function
                    
                    [dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT(l,1:3), dataGroupPool(subj,Session).hrbyStimbintACSLagMain_RESNORM(l,1),dataGroupPool(subj,Session).hrbyStimbintACSLagMain_RESIDUAL(l,1:18),dataGroupPool(subj,Session).hrbyStimbintACSLagMain_prefPhase(l,1)] = fitCos2RatebyPhase ([dataGroupPool(subj,Session).hrbyStimbintACSLagMain(l,:) dataGroupPool(subj,Session).hrbyStimbintACSLagMain(l,:)]',X,params);
                    dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT_yAll(l,1:9) = dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT(l,2) + dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT(l,3).*(cos(X + dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT(l,1))); % if you want to plot the predicted function
                catch
                    errorPL = [subj Session l];
                dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT(l,:) = nan;
                dataGroupPool(subj,Session).hrbytACSbinStimLagMain_RESNORM(l,:)= nan;
                dataGroupPool(subj,Session).hrbytACSbinStimLagMain_RESIDUAL(l,:)= nan;
                dataGroupPool(subj,Session).hrbytACSbinStimLagMain_prefPhase(l,:) = nan;
                dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT_yAll(l,:) = nan;
                dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT(l,:)= nan;
                dataGroupPool(subj,Session).hrbyStimbintACSLagMain_RESNORM(l,:)= nan;
                dataGroupPool(subj,Session).hrbyStimbintACSLagMain_RESIDUAL(l,:)= nan;
                dataGroupPool(subj,Session).hrbyStimbintACSLagMain_prefPhase(l,:)= nan;
                dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT_yAll(l,:) = nan;

                end
            end
            %fit only by gap phas
            [dataGroupPool(subj,Session).hrbyStimbinMain_FIT, dataGroupPool(subj,Session).hrbyStimbinMain_RESNORM,dataGroupPool(subj,Session).hrbyStimbinMain_RESIDUAL,dataGroupPool(subj,Session).hrbyStimbinMain_prefPhase] = fitCos2RatebyPhase ([dataGroupPool(subj,Session).hrbyStimbinMain dataGroupPool(subj,Session).hrbyStimbinMain]',X,params);
            dataGroupPool(subj,Session).hrbyStimbinMain_FIT_yAll = dataGroupPool(subj,Session).hrbyStimbinMain_FIT(2) + dataGroupPool(subj,Session).hrbyStimbinMain_FIT(3).*(cos(X + dataGroupPool(subj,Session).hrbyStimbinMain_FIT(1))); % if you want to plot the predicted function
            
            [dataGroupPool(subj,Session).hrbyStimbinsham_FIT, dataGroupPool(subj,Session).hrbyStimbinsham_RESNORM,dataGroupPool(subj,Session).hrbyStimbinsham_RESIDUAL,dataGroupPool(subj,Session).hrbyStimbinsham_prefPhase] = fitCos2RatebyPhase ([dataGroupPool(subj,Session).hrbyStimbinsham dataGroupPool(subj,Session).hrbyStimbinsham]',X,params);
            dataGroupPool(subj,Session).hrbyStimbinsham_FIT_yAll = dataGroupPool(subj,Session).hrbyStimbinsham_FIT(2) + dataGroupPool(subj,Session).hrbyStimbinsham_FIT(3).*(cos(X + dataGroupPool(subj,Session).hrbyStimbinsham_FIT(1))); % if you want to plot the predicted function
            
            [dataGroupPool(subj,Session).hrbytACSbinMain_FIT, dataGroupPool(subj,Session).hrbytACSbinMain_RESNORM,dataGroupPool(subj,Session).hrbytACSbinMain_RESIDUAL,dataGroupPool(subj,Session).hrbytACSbinMain_prefPhase] = fitCos2RatebyPhase ([dataGroupPool(subj,Session).hrbytACSbinMain dataGroupPool(subj,Session).hrbytACSbinMain]',X,params);
            dataGroupPool(subj,Session).hrbytACSbinMain_FIT_yAll = dataGroupPool(subj,Session).hrbytACSbinMain_FIT(2) + dataGroupPool(subj,Session).hrbytACSbinMain_FIT(3).*(cos(X + dataGroupPool(subj,Session).hrbytACSbinMain_FIT(1))); % if you want to plot the predicted function
            
            [dataGroupPool(subj,Session).hrbytACSbinsham_FIT, dataGroupPool(subj,Session).hrbytACSbinsham_RESNORM,dataGroupPool(subj,Session).hrbytACSbinsham_RESIDUAL,dataGroupPool(subj,Session).hrbytACSbinsham_prefPhase] = fitCos2RatebyPhase ([dataGroupPool(subj,Session).hrbytACSbinsham dataGroupPool(subj,Session).hrbytACSbinsham]',X,params);
            dataGroupPool(subj,Session).hrbytACSbinsham_FIT_yAll = dataGroupPool(subj,Session).hrbytACSbinsham_FIT(2) + dataGroupPool(subj,Session).hrbytACSbinsham_FIT(3).*(cos(X + dataGroupPool(subj,Session).hrbytACSbinsham_FIT(1))); % if you want to plot the predicted function
            
            dataGroupPool(subj,Session).HRall = [dataGroupPool(subj,Session).hrbyStimbinsham_FIT(2) dataGroupPool(subj,Session).hrbyStimbinMain_FIT(:,2) dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT(:,2)'];
            dataGroupPool(subj,Session).HRalltacs = [dataGroupPool(subj,Session).hrbytACSbinsham_FIT(2) dataGroupPool(subj,Session).hrbytACSbinMain_FIT(:,2) dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT(:,2)'];
            
            dataGroupPool(subj,Session).Ampall = [dataGroupPool(subj,Session).hrbyStimbinsham_FIT(3) dataGroupPool(subj,Session).hrbyStimbinMain_FIT(:,3) dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT(:,3)'];
            dataGroupPool(subj,Session).Ampalltacs = [dataGroupPool(subj,Session).hrbytACSbinsham_FIT(3) dataGroupPool(subj,Session).hrbytACSbinMain_FIT(:,3) dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT(:,3)'];
            close all
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %create surrogate datasets by shuffeling the tACS phase lag
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            gapAUDtACSlagMainS = nan(length(gapAUDtACSlagMain),NrIte);
            %               for ite =1:NrIte
            ite=0;
            while ite<NrIte+1
                try
                    ite=ite+1;
                    disp(ite)
                    neworder = randperm(length(gapAUDtACSlagMain));
                    gapAUDtACSlagMainS(:,ite) = gapAUDtACSlagMain(neworder);
                    %gapAUDtACSlagMain = wrapTo2Pi(gapAUD_tACSlag(gapStimcond==1)); %phase lag between tACS and audio signal
                    for gapBin = 1:length(phasebinGap)-1
                        currGapbinStim = find(gapAUDphaseMain>=phasebinGap(gapBin)&gapAUDphaseMain<phasebinGap(gapBin+1))'; %gaps im current AUD bin
                        length(currGapbinStim)
                        for b=1:length(lagBins)-1 %separate by tACS-AUD phase lag
                            if lagBins1(b)>lagBins2(b)
                                currGapStimbinSlag = find((gapAUDphaseMain>=phasebinGap(gapBin)&gapAUDphaseMain<phasebinGap(gapBin+1))&...
                                    (gapAUDtACSlagMainS(:,ite)>=lagBins1(b) | gapAUDtACSlagMainS(:,ite)<lagBins2(b)))';
                            else
                                currGapStimbinSlag =find((gapAUDphaseMain>=phasebinGap(gapBin)&gapAUDphaseMain<phasebinGap(gapBin+1))&...
                                    (gapAUDtACSlagMainS(:,ite)>=lagBins1(b)& gapAUDtACSlagMainS(:,ite)<lagBins2(b)))';
                            end
                            length(currGapStimbinSlag)
                            dataGroupPool(subj,Session).hrbyStimbintACSLagMainS(ite,b,gapBin) = mean(trCurrGapTrackMain(currGapStimbinSlag));
                            
                                                    figure,
                                                    subplot(1,2,1), circ_plot(gapAUDphaseMain(currGapStimbinSlag))
                                                    title(['AUDIO gap bin: ' num2str(gapBin) ' lag bin' num2str(b)])
                                                    subplot(1,2,2), circ_plot(gapAUDtACSlagMainS(currGapStimbinSlag,ite))
                                                    title('lag')
                        end
                    end
                    for l = 1:length(lagBins)-1
                        %                       try
                        %fit by gap phase and stim_tacs phase lag
                        [dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FITS(ite,l,:), dataGroupPool(subj,Session).hrbyStimbintACSLagMain_RESNORMS(ite,l,:),dataGroupPool(subj,Session).hrbyStimbintACSLagMain_RESIDUALS(ite,l,:),dataGroupPool(subj,Session).hrbyStimbintACSLagMain_prefPhaseS(ite,l,:)] = fitCos2RatebyPhase ([squeeze(dataGroupPool(subj,Session).hrbyStimbintACSLagMainS(ite,l,:)); squeeze(dataGroupPool(subj,Session).hrbyStimbintACSLagMainS(ite,l,:))],X,params);
                        dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT_yAllS(ite,l,:) = dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FITS(ite,l,2) + dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FITS(ite,l,3).*(cos(X + dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FITS(ite,l,1))); % if you want to plot the predicted function
                        %                       catch
                        %                       end
                    end
                    
                    %shuffle hit rates now and fit only by audio phase
                    trCurrGapTrackshamS = trCurrGapTracksham(randperm(length(trCurrGapTracksham)));
                    for gapBin = 1:length(phasebinGap)-1
                        currGapbinStimshamS = find(gapAUDphasesham>=phasebinGap(gapBin)&gapAUDphasesham<phasebinGap(gapBin+1))';
                        %    length(currGapbinStimshamS)
                        dataGroupPool(subj,Session).hrbyStimbinshamS(ite,gapBin) = mean(trCurrGapTrackshamS(currGapbinStimshamS));
                    end
                    %fit only by gap phas
                    [dataGroupPool(subj,Session).hrbyStimbinsham_FITS(ite,:), dataGroupPool(subj,Session).hrbyStimbinsham_RESNORMS(ite,:),dataGroupPool(subj,Session).hrbyStimbinsham_RESIDUALS(ite,:),dataGroupPool(subj,Session).hrbyStimbinsham_prefPhaseS(ite,:)] = fitCos2RatebyPhase ([dataGroupPool(subj,Session).hrbyStimbinshamS(ite,:) dataGroupPool(subj,Session).hrbyStimbinshamS(ite,:)]',X,params);
                    dataGroupPool(subj,Session).hrbyStimbinsham_FIT_yAllS(ite,:) = dataGroupPool(subj,Session).hrbyStimbinsham_FITS(ite,2) + dataGroupPool(subj,Session).hrbyStimbinsham_FITS(ite,3).*(cos(X + dataGroupPool(subj,Session).hrbyStimbinsham_FITS(ite,1))); % if you want to plot the predicted function
                    close all
                catch
                    ite =ite-1;
                end
            end
        
        end
    end
end

%%%%%%%%%%%SAME ANALYSIS AS BEFORE BUT POOLING DATA ACROSS SESSIONS
for subj=1:length(SUBJlist)%[1 3:6 8:26 28:38 40:length(SUBJlist)]%loop across subjects
    close all
    currSUBJ = SUBJlist{subj,1}; % current subject
    fprintf('%s%s\n','... working on subject ', currSUBJ)
    
    NrSess   = SUBJlist{subj,2};
    gapDetect = [];
    gaptACSangle = [];
    gapAUDangle  = [];
    gapAUD_tACSlag = [];
    gapStimcond = [];
    
    if ~isempty(dataGroupallGaps(subj,1).gapDetect)&&~isempty(dataGroupallGaps(subj,2).gapDetect) %check that there is data for both sessions
        %concatenate data from 2 sessions
        gapDetect = [dataGroupallGaps(subj,1).gapDetect; dataGroupallGaps(subj,2).gapDetect];
        gaptACSangle = [dataGroupallGaps(subj,1).gaptACSangle dataGroupallGaps(subj,2).gaptACSangle]';
        gapAUDangle = [dataGroupallGaps(subj,1).gapAUDangle dataGroupallGaps(subj,2).gapAUDangle]';
        gapAUD_tACSlag = [dataGroupallGaps(subj,1).gapAUD_tACSlag dataGroupallGaps(subj,2).gapAUD_tACSlag]';
        gapStimcond = [dataGroupallGaps(subj,1).gapStimcond; dataGroupallGaps(subj,2).gapStimcond];
        
        %first separate sham from real
        %main
        tACS1angleGapMain  = wrapTo2Pi(gaptACSangle(gapStimcond==1)); %tACS phase at gap onset
        trCurrGapTrackMain = gapDetect(gapStimcond==1); %gap info from eeg triggers
        gapAUDtACSlagMain  = wrapTo2Pi(gapAUD_tACSlag(gapStimcond==1)); %phase lag between tACS and audio signal
        gapAUDphaseMain    = gapAUDangle(gapStimcond==1); %FM phase at gap onset, 1=real, 2=sham
        %sham
        tACS1angleGapsham  = wrapTo2Pi(gaptACSangle(gapStimcond==2)); %tACS phase at gap onset
        trCurrGapTracksham = gapDetect(gapStimcond==2); %gap info from eeg triggers
        gapAUDtACSlagsham  = wrapTo2Pi(gapAUD_tACSlag(gapStimcond==2)); %phase lag between tACS and audio signal
        gapAUDphasesham    = gapAUDangle(gapStimcond==2); %FM phase at gap onset
        
        for gapBin = 1:length(phasebinGap)-1
            currGapbintACS  = find(tACS1angleGapMain>=phasebinGap(gapBin)&tACS1angleGapMain<phasebinGap(gapBin+1))'; %gaps in current tACS phase bin
            currGapbinStim  = find(gapAUDphaseMain>=phasebinGap(gapBin)&gapAUDphaseMain<phasebinGap(gapBin+1))'; %gaps im current AUD bin
            length(currGapbintACS);
            length(currGapbinStim);
            for b=1:length(lagBins)-1 %separate by tACS-AUD phase lag
                if lagBins1(b)>lagBins2(b)
                    currGaptACSbinSlag = find((tACS1angleGapMain>=phasebinGap(gapBin)&tACS1angleGapMain<phasebinGap(gapBin+1))&...
                        (gapAUDtACSlagMain>=lagBins1(b) | gapAUDtACSlagMain<lagBins2(b)))';
                    currGapStimbinSlag = find((gapAUDphaseMain>=phasebinGap(gapBin)&gapAUDphaseMain<phasebinGap(gapBin+1))&...
                        (gapAUDtACSlagMain>=lagBins1(b) | gapAUDtACSlagMain<lagBins2(b)))';
                else
                    currGaptACSbinSlag = find((tACS1angleGapMain>=phasebinGap(gapBin)&tACS1angleGapMain<phasebinGap(gapBin+1))&...
                        (gapAUDtACSlagMain>=lagBins1(b)&gapAUDtACSlagMain<lagBins2(b)))';
                    currGapStimbinSlag =find((gapAUDphaseMain>=phasebinGap(gapBin)&gapAUDphaseMain<phasebinGap(gapBin+1))&...
                        (gapAUDtACSlagMain>=lagBins1(b)& gapAUDtACSlagMain<lagBins2(b)))';
                end
                length(currGaptACSbinSlag);
                length(currGapStimbinSlag);
                dataGroupPool(subj,3).hrbytACSbinStimLagMainALL(b,gapBin) = mean(trCurrGapTrackMain(currGaptACSbinSlag));
                dataGroupPool(subj,3).hrbyStimbintACSLagMainALL(b,gapBin) = mean(trCurrGapTrackMain(currGapStimbinSlag));
                %
                %                                         figure, subplot(2,2,1), circ_plot(tACS1angleGapMain(currGaptACSbinSlag))
                %                                         title(['TACS gap bin: ' num2str(gapBin) ' lag bin' num2str(b)])
                %                                         subplot(2,2,2), circ_plot(gapAUDtACSlagMain(currGaptACSbinSlag))
                %                                         title('lag')
                %                                         subplot(2,2,3), circ_plot(gapAUDphaseMain(currGapStimbinSlag))
                %                                         title(['AUDIO gap bin: ' num2str(gapBin) ' lag bin' num2str(b)])
                %                                         subplot(2,2,4), circ_plot(gapAUDtACSlagMain(currGapStimbinSlag))
                %                                         title('lag')
            end
            %separate by both FM and tACS signal
            for gapBin2 = 1:length(phasebinGap)-1
                currGapbinAT = find((tACS1angleGapMain>=phasebinGap(gapBin)&tACS1angleGapMain<phasebinGap(gapBin+1))&...
                    (gapAUDphaseMain>=phasebinGap(gapBin2)&gapAUDphaseMain<phasebinGap(gapBin2+1)))';
                
                %  length(currGapbinAT)
                dataGroupPool(subj,3).hrbyStimbytACSbinMainALL(gapBin2,gapBin) = mean(trCurrGapTrackMain(currGapbinAT));
                dataGroupPool(subj,3).hrbyStimbytACSbinMainNrALL(gapBin2,gapBin) = length(currGapbinAT);
            end
            dataGroupPool(subj,3).hrbyStimbinMainALL(gapBin) = mean(trCurrGapTrackMain(currGapbinStim));
            dataGroupPool(subj,3).hrbytACSbinMainALL(gapBin) = mean(trCurrGapTrackMain(currGapbintACS));
        end
        
        for gapBin = 1:length(phasebinGap)-1
            
            currGapbintACSsham = find(tACS1angleGapsham>=phasebinGap(gapBin)&tACS1angleGapsham<phasebinGap(gapBin+1))';
            currGapbinStimsham = find(gapAUDphasesham>=phasebinGap(gapBin)&gapAUDphasesham<phasebinGap(gapBin+1))';
            %                 length(currGapbintACSsham)
            %                 length(currGapbinStimsham)
            dataGroupPool(subj,3).hrbyStimbinshamALL(gapBin) = mean(trCurrGapTracksham(currGapbinStimsham));
            dataGroupPool(subj,3).hrbytACSbinshamALL(gapBin) = mean(trCurrGapTracksham(currGapbintACSsham));
            
            %separate by both FM and tACS signal
            for gapBin2 = 1:length(phasebinGap)-1
                currGapbinATs = find((tACS1angleGapsham>=phasebinGap(gapBin)&tACS1angleGapsham<phasebinGap(gapBin+1))&...
                    (gapAUDphasesham>=phasebinGap(gapBin2)&gapAUDphasesham<phasebinGap(gapBin2+1)))';
                %  length(currGapbinATs)
                dataGroupPool(subj,3).hrbyStimbytACSbinshamALL(gapBin2,gapBin) = mean(trCurrGapTracksham(currGapbinATs));
                dataGroupPool(subj,3).hrbyStimbytACSbinshamNrALL(gapBin2,gapBin) = length(currGapbinATs);
                
            end
        end
        
%         %fit Cosine function
%         X         = 0:2*pi/NrGapBins:2*pi-2*pi/NrGapBins; %this is my phase vector
%         X         = X+2*pi/NrGapBins/2;
%         lag       = pi;
%         intercept = 0.5;
%         amp       = 0.5;
%         params    = [lag intercept amp];
        for l = 1:length(lagBins)-1
            try
                %fit by gap phase and stim_tacs phase lag
                [dataGroupPool(subj,3).hrbytACSbinStimLagMain_FITALL(l,:), dataGroupPool(subj,3).hrbytACSbinStimLagMain_RESNORMALL(l,:),dataGroupPool(subj,3).hrbytACSbinStimLagMain_RESIDUALALL(l,:),dataGroupPool(subj,3).hrbytACSbinStimLagMain_prefPhaseALL(l,:)] = fitCos2RatebyPhase ([dataGroupPool(subj,3).hrbytACSbinStimLagMainALL(l,:) dataGroupPool(subj,3).hrbytACSbinStimLagMainALL(l,:)]',X,params);
                dataGroupPool(subj,3).hrbytACSbinStimLagMain_FIT_yAllALL(l,:) = dataGroupPool(subj,3).hrbytACSbinStimLagMain_FITALL(l,2) + dataGroupPool(subj,3).hrbytACSbinStimLagMain_FITALL(l,3).*(cos(X + dataGroupPool(subj,3).hrbytACSbinStimLagMain_FITALL(l,1))); % if you want to plot the predicted function
                
                [dataGroupPool(subj,3).hrbyStimbintACSLagMain_FITALL(l,:), dataGroupPool(subj,3).hrbyStimbintACSLagMain_RESNORMALL(l,:),dataGroupPool(subj,3).hrbyStimbintACSLagMain_RESIDUALALL(l,:),dataGroupPool(subj,3).hrbyStimbintACSLagMain_prefPhaseALL(l,:)] = fitCos2RatebyPhase ([dataGroupPool(subj,3).hrbyStimbintACSLagMainALL(l,:) dataGroupPool(subj,3).hrbyStimbintACSLagMainALL(l,:)]',X,params);
                dataGroupPool(subj,3).hrbyStimbintACSLagMain_FIT_yAllALL(l,:) = dataGroupPool(subj,3).hrbyStimbintACSLagMain_FITALL(l,2) + dataGroupPool(subj,3).hrbyStimbintACSLagMain_FITALL(l,3).*(cos(X + dataGroupPool(subj,3).hrbyStimbintACSLagMain_FITALL(l,1))); % if you want to plot the predicted function
            catch
                errorPL(subj,l) = 1;
                dataGroupPool(subj,3).hrbytACSbinStimLagMain_FITALL(l,:) = nan;
                dataGroupPool(subj,3).hrbytACSbinStimLagMain_RESNORMALL(l,:)= nan;
                dataGroupPool(subj,3).hrbytACSbinStimLagMain_RESIDUALALL(l,:)= nan;
                dataGroupPool(subj,3).hrbytACSbinStimLagMain_prefPhaseALL(l,:) = nan;
                dataGroupPool(subj,3).hrbytACSbinStimLagMain_FIT_yAllALL(l,:) = nan;
                dataGroupPool(subj,3).hrbyStimbintACSLagMain_FITALL(l,:)= nan;
                dataGroupPool(subj,3).hrbyStimbintACSLagMain_RESNORMALL(l,:)= nan;
                dataGroupPool(subj,3).hrbyStimbintACSLagMain_RESIDUALALL(l,:)= nan;
                dataGroupPool(subj,3).hrbyStimbintACSLagMain_prefPhaseALL(l,:)= nan;
                dataGroupPool(subj,3).hrbyStimbintACSLagMain_FIT_yAllALL(l,:) = nan;
            end
        end
        %fit only by gap phas
        [dataGroupPool(subj,3).hrbyStimbinMain_FITALL, dataGroupPool(subj,3).hrbyStimbinMain_RESNORMALL,dataGroupPool(subj,3).hrbyStimbinMain_RESIDUALALL,dataGroupPool(subj,3).hrbyStimbinMain_prefPhaseALL] = fitCos2RatebyPhase ([dataGroupPool(subj,3).hrbyStimbinMainALL dataGroupPool(subj,3).hrbyStimbinMainALL]',X,params);
        dataGroupPool(subj,3).hrbyStimbinMain_FIT_yAllALL = dataGroupPool(subj,3).hrbyStimbinMain_FITALL(2) + dataGroupPool(subj,3).hrbyStimbinMain_FITALL(3).*(cos(X + dataGroupPool(subj,3).hrbyStimbinMain_FITALL(1))); % if you want to plot the predicted function
        
        [dataGroupPool(subj,3).hrbyStimbinsham_FITALL, dataGroupPool(subj,3).hrbyStimbinsham_RESNORMALL,dataGroupPool(subj,3).hrbyStimbinsham_RESIDUALALL,dataGroupPool(subj,3).hrbyStimbinsham_prefPhaseALL] = fitCos2RatebyPhase ([dataGroupPool(subj,3).hrbyStimbinshamALL dataGroupPool(subj,3).hrbyStimbinshamALL]',X,params);
        dataGroupPool(subj,3).hrbyStimbinsham_FIT_yAllALL = dataGroupPool(subj,3).hrbyStimbinsham_FITALL(2) + dataGroupPool(subj,3).hrbyStimbinsham_FITALL(3).*(cos(X + dataGroupPool(subj,3).hrbyStimbinsham_FITALL(1))); % if you want to plot the predicted function
        
        [dataGroupPool(subj,3).hrbytACSbinMain_FITALL, dataGroupPool(subj,3).hrbytACSbinMain_RESNORMALL,dataGroupPool(subj,3).hrbytACSbinMain_RESIDUALALL,dataGroupPool(subj,3).hrbytACSbinMain_prefPhaseALL] = fitCos2RatebyPhase ([dataGroupPool(subj,3).hrbytACSbinMainALL dataGroupPool(subj,3).hrbytACSbinMainALL]',X,params);
        dataGroupPool(subj,3).hrbytACSbinMain_FIT_yAllALL = dataGroupPool(subj,3).hrbytACSbinMain_FITALL(2) + dataGroupPool(subj,3).hrbytACSbinMain_FITALL(3).*(cos(X + dataGroupPool(subj,3).hrbytACSbinMain_FITALL(1))); % if you want to plot the predicted function
        
        [dataGroupPool(subj,3).hrbytACSbinsham_FITALL, dataGroupPool(subj,3).hrbytACSbinsham_RESNORMALL,dataGroupPool(subj,3).hrbytACSbinsham_RESIDUALALL,dataGroupPool(subj,3).hrbytACSbinsham_prefPhaseALL] = fitCos2RatebyPhase ([dataGroupPool(subj,3).hrbytACSbinshamALL dataGroupPool(subj,3).hrbytACSbinshamALL]',X,params);
        dataGroupPool(subj,3).hrbytACSbinsham_FIT_yAllALL = dataGroupPool(subj,3).hrbytACSbinsham_FITALL(2) + dataGroupPool(subj,3).hrbytACSbinsham_FITALL(3).*(cos(X + dataGroupPool(subj,3).hrbytACSbinsham_FITALL(1))); % if you want to plot the predicted function
        
        dataGroupPool(subj,3).HRallALL = [dataGroupPool(subj,3).hrbyStimbinsham_FITALL(2) dataGroupPool(subj,3).hrbyStimbinMain_FITALL(:,2) dataGroupPool(subj,3).hrbyStimbintACSLagMain_FITALL(:,2)'];
        dataGroupPool(subj,3).HRalltacsALL = [dataGroupPool(subj,3).hrbytACSbinsham_FITALL(2) dataGroupPool(subj,3).hrbytACSbinMain_FITALL(:,2) dataGroupPool(subj,3).hrbytACSbinStimLagMain_FITALL(:,2)'];
        
        dataGroupPool(subj,3).AmpallALL = [dataGroupPool(subj,3).hrbyStimbinsham_FITALL(3) dataGroupPool(subj,3).hrbyStimbinMain_FITALL(:,3) dataGroupPool(subj,3).hrbyStimbintACSLagMain_FITALL(:,3)'];
        dataGroupPool(subj,3).AmpalltacsALL = [dataGroupPool(subj,3).hrbytACSbinsham_FITALL(3) dataGroupPool(subj,3).hrbytACSbinMain_FITALL(:,3) dataGroupPool(subj,3).hrbytACSbinStimLagMain_FITALL(:,3)'];
        close all
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %create surrogate datasets by shuffeling the tACS phase lag
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        gapAUDtACSlagMainS = nan(length(gapAUDtACSlagMain),NrIte);
        %               for ite =1:NrIte
        ite=0;
        while ite<NrIte+1
            try
                ite=ite+1;
                disp(ite)
                neworder = randperm(length(gapAUDtACSlagMain));
                gapAUDtACSlagMainS(:,ite) = gapAUDtACSlagMain(neworder);
                %gapAUDtACSlagMain = wrapTo2Pi(gapAUD_tACSlag(gapStimcond==1)); %phase lag between tACS and audio signal
                for gapBin = 1:length(phasebinGap)-1
                    currGapbinStim = find(gapAUDphaseMain>=phasebinGap(gapBin)&gapAUDphaseMain<phasebinGap(gapBin+1))'; %gaps im current AUD bin
                    length(currGapbinStim);
                    for b=1:length(lagBins)-1 %separate by tACS-AUD phase lag
                        if lagBins1(b)>lagBins2(b)
                            currGapStimbinSlag = find((gapAUDphaseMain>=phasebinGap(gapBin)&gapAUDphaseMain<phasebinGap(gapBin+1))&...
                                (gapAUDtACSlagMainS(:,ite)>=lagBins1(b) | gapAUDtACSlagMainS(:,ite)<lagBins2(b)))';
                        else
                            currGapStimbinSlag =find((gapAUDphaseMain>=phasebinGap(gapBin)&gapAUDphaseMain<phasebinGap(gapBin+1))&...
                                (gapAUDtACSlagMainS(:,ite)>=lagBins1(b)& gapAUDtACSlagMainS(:,ite)<lagBins2(b)))';
                        end
                        length(currGapStimbinSlag);
                        dataGroupPool(subj,3).hrbyStimbintACSLagMainSALL(ite,b,gapBin) = mean(trCurrGapTrackMain(currGapStimbinSlag));
                        
                        %                         figure,
                        %                         subplot(1,2,1), circ_plot(gapAUDphaseMain(currGapStimbinSlag))
                        %                         title(['AUDIO gap bin: ' num2str(gapBin) ' lag bin' num2str(b)])
                        %                         subplot(1,2,2), circ_plot(gapAUDtACSlagMainS(currGapStimbinSlag,ite))
                        %                         title('lag')
                    end
                end
                for l = 1:length(lagBins)-1
                    %                       try
                    %fit by gap phase and stim_tacs phase lag
                    [dataGroupPool(subj,3).hrbyStimbintACSLagMain_FITSALL(ite,l,:), dataGroupPool(subj,3).hrbyStimbintACSLagMain_RESNORMSALL(ite,l,:),dataGroupPool(subj,3).hrbyStimbintACSLagMain_RESIDUALSALL(ite,l,:),dataGroupPool(subj,3).hrbyStimbintACSLagMain_prefPhaseSALL(ite,l,:)] = fitCos2RatebyPhase ([squeeze(dataGroupPool(subj,3).hrbyStimbintACSLagMainSALL(ite,l,:)); squeeze(dataGroupPool(subj,3).hrbyStimbintACSLagMainSALL(ite,l,:))],X,params);
                    dataGroupPool(subj,3).hrbyStimbintACSLagMain_FIT_yAllSALL(ite,l,:) = dataGroupPool(subj,3).hrbyStimbintACSLagMain_FITSALL(ite,l,2) + dataGroupPool(subj,3).hrbyStimbintACSLagMain_FITSALL(ite,l,3).*(cos(X + dataGroupPool(subj,3).hrbyStimbintACSLagMain_FITSALL(ite,l,1))); % if you want to plot the predicted function
                    %                       catch
                    %                       end
                end
                
                %shuffle hit rates now and fit only by audio phase
                trCurrGapTrackshamS = trCurrGapTracksham(randperm(length(trCurrGapTracksham)));
                for gapBin = 1:length(phasebinGap)-1
                    currGapbinStimshamS = find(gapAUDphasesham>=phasebinGap(gapBin)&gapAUDphasesham<phasebinGap(gapBin+1))';
                    %    length(currGapbinStimshamS)
                    dataGroupPool(subj,3).hrbyStimbinshamSALL(ite,gapBin) = mean(trCurrGapTrackshamS(currGapbinStimshamS));
                end
                %fit only by gap phas
                [dataGroupPool(subj,3).hrbyStimbinsham_FITSALL(ite,:), dataGroupPool(subj,3).hrbyStimbinsham_RESNORMSALL(ite,:),dataGroupPool(subj,3).hrbyStimbinsham_RESIDUALSALL(ite,:),dataGroupPool(subj,3).hrbyStimbinsham_prefPhaseSALL(ite,:)] = fitCos2RatebyPhase ([dataGroupPool(subj,3).hrbyStimbinshamSALL(ite,:) dataGroupPool(subj,3).hrbyStimbinshamSALL(ite,:)]',X,params);
                dataGroupPool(subj,3).hrbyStimbinsham_FIT_yAllSALL(ite,:) = dataGroupPool(subj,3).hrbyStimbinsham_FITSALL(ite,2) + dataGroupPool(subj,3).hrbyStimbinsham_FITSALL(ite,3).*(cos(X + dataGroupPool(subj,3).hrbyStimbinsham_FITSALL(ite,1))); % if you want to plot the predicted function
                close all
            catch
                ite =ite-1;
            end
        end
        %               %new strucut with multiple bins and overlaps
        %               dataLast(subj,Session).data = dataGroupPool(subj,Session);
    end
end

save (fullfile(OUTdir,'dataGroupPooltACS_eLifeRev2.mat'), 'dataGroupPool', '-v7.3')
load(fullfile(OUTdir,'dataGroupPooltACS_eLifeRev2.mat'),'dataGroupPool')


%fit cos to tacs lag for session 2 as well
tee=1;
figure
for Session=1:2
    
    for s = 1:size(dataGroupPool,1)
        try
            groupmatrixFin.data(s,:,Session) = dataGroupPool(s,Session).hrbyStimbintACSLagMain_FIT(:,3)';
            X         = 0:2*pi/NrAudtACSlag:2*pi-2*pi/NrAudtACSlag; %this is my phase vector
            Xplot     = (0:2*pi/1000:2*pi-2*pi/NrAudtACSlag)+2*pi/NrAudtACSlag/2; %this is my phase vector
            X         = X+2*pi/NrAudtACSlag/2;
            lag       = pi;
            intercept = 0.25;
            amp       = 0.1;
            params    = [lag intercept amp];
            
            TEMPdata = groupmatrixFin.data(s,:,Session);
            [groupmatrixFin.dataFIT(s,:,Session), ~,~,groupmatrixFin.dataprefPhase(s,Session)] = fitCos2RatebyPhase ([TEMPdata TEMPdata]',X,params);
            groupmatrixFin.AlignCos_FIT_yAll(s,Session,:) = groupmatrixFin.dataFIT(s,2,Session) + groupmatrixFin.dataFIT(s,3,Session).*(cos(Xplot + groupmatrixFin.dataFIT(s,1,Session)))'; % if you want to plot the predicted function
            
            
            %run the fit for the permuted datasets
            for ite =1:1000
                groupmatrixFin.dataS(s,ite,:,Session)=dataGroupPool(s,Session).hrbyStimbintACSLagMain_FITS(ite,:,3);
                TEMPdataS = squeeze(groupmatrixFin.dataS(s,ite,:,Session))';
                [groupmatrixFin.dataFITS(s,ite,:,Session), ~,~,groupmatrixFin.dataprefPhaseS(s,ite,Session)] = fitCos2RatebyPhase ([TEMPdataS TEMPdataS]',X,params);
                groupmatrixFin.AlignCos_FIT_yAllS(s,ite,:,Session) = groupmatrixFin.dataFITS(s,ite,2,Session) + groupmatrixFin.dataFITS(s,ite,3,Session).*(cos(X + groupmatrixFin.dataFITS(s,ite,1,Session))); % if you want to plot the predicted function
            end
            subplot(7,6,s)
            plot(X,squeeze(groupmatrixFin.data(s,:,Session)))
            hold on
            plot(Xplot,squeeze(groupmatrixFin.AlignCos_FIT_yAll(s,Session,:)))
            
            title(['Subj' num2str(s) '_S' num2str(Session)])
        catch
            probsu = s;
            tee=tee+1;
            groupmatrixFin.dataFIT(s,:,Session)=nan;
            groupmatrixFin.dataprefPhase(s,Session)=nan;
        end
        
    end
end
%only plot
figure
for Session=1:2
    
    for s = 1:size(dataGroupPool,1)
        try
            subplot(7,6,s)
            hold on
            plot(X,squeeze(groupmatrixFin.data(s,:,Session)))
            hold on
            plot(Xplot,squeeze(groupmatrixFin.AlignCos_FIT_yAll(s,Session,:)))
            
            title(['Subj' num2str(s) '_S' num2str(Session)])
        catch
           
        end
        
    end
end
%temporal save in case it crashes
save eLifeRev2

for Session=1:2
    figure
    for s = 1:size(dataGroupPool,1)
            subplot(7,6,s)
            plot([dataGroupPool(s,Session).hrbyStimbintACSLagMain dataGroupPool(s,Session).hrbyStimbintACSLagMain]')
            hold on
             plot([dataGroupPool(s,Session).hrbyStimbinsham dataGroupPool(s,Session).hrbyStimbinsham]')
             plot([dataGroupPool(s,Session).hrbyStimbintACSLagMain_FIT_yAll dataGroupPool(s,Session).hrbyStimbintACSLagMain_FIT_yAll]','--')
             plot([dataGroupPool(s,Session).hrbyStimbinsham_FIT_yAll dataGroupPool(s,Session).hrbyStimbinsham_FIT_yAll]','--')    
            title(['Subj' num2str(s) '_S' num2str(Session)])
            ylim([0 1])
    if s ==1
        legend
    end
        
    end
end

dataAmpS       = nan(length(SUBJlist),NrIte,2);
dataAmpOrig    = nan(length(SUBJlist),2);
dataPhaseOrig  = nan(length(SUBJlist),2);
dataAmpOrigInt = nan(length(SUBJlist),2);
dataPhaseOrigS  = nan(length(SUBJlist),NrIte,2);

for Session =1:2
    dataAmpS(:,:,Session)     = squeeze(groupmatrixFin.dataFITS(:,:,3,Session));
    dataAmpOrig(:,Session)    = squeeze(groupmatrixFin.dataFIT(:,3,Session));
    dataPhaseOrig(:,Session)  = squeeze(groupmatrixFin.dataprefPhase(:,Session));
    dataPhaseOrigS(:,:,Session)  = squeeze(groupmatrixFin.dataprefPhaseS(:,:,Session));
    dataAmpOrigInt(:,Session) = squeeze(groupmatrixFin.dataFIT(:,2,Session));  
end

%compute subject z-score relative to permuted datasets
dataAmpZscoreOrig = (dataAmpOrig-squeeze(mean(dataAmpS,2)))./squeeze(std(dataAmpS,[],2));
%FIT THE COSINE POOLING DATA ACROSS SESSIONS
sp2 = 1;
%plot amp by bin and overlap
tee=1;
figure
sp=1;

save eLifeRev2

for s = 1:size(dataGroupPool,1)
    try
        X         = 0:2*pi/NrAudtACSlag:2*pi-2*pi/NrAudtACSlag; %this is my phase vector
        X         = X+2*pi/NrAudtACSlag/2;
        lag       = pi;
        intercept = 0.25;
        amp       = 0.1;
        params    = [lag intercept amp];
        
        TEMPdata = dataGroupPool(s,3).hrbyStimbintACSLagMain_FITALL(:,3)';
        [groupmatrixFin.dataFITALL(s,:), ~,~,groupmatrixFin.dataprefPhaseAll(s)] = fitCos2RatebyPhase ([TEMPdata TEMPdata]',X,params);
        groupmatrixFin.AlignCos_FIT_yAllALL(s,:) = groupmatrixFin.dataFITALL(s,2) + groupmatrixFin.dataFITALL(s,3).*(cos(X + groupmatrixFin.dataFITALL(s,1))); % if you want to plot the predicted function
        groupmatrixFin.dataALL(s,:) = TEMPdata;
        
        %run the fit for the permuted datasets
        for ite =1:1000
            TEMPdataS = dataGroupPool(s,3).hrbyStimbintACSLagMain_FITSALL(ite,:,3);
            groupmatrixFin.dataSALL(s,ite,:) = TEMPdataS;
            
            [groupmatrixFin.dataFITSALL(s,ite,:), ~,~,groupmatrixFin.dataprefPhaseSALL(s,ite)] = fitCos2RatebyPhase ([TEMPdataS TEMPdataS]',X,params);
            groupmatrixFin.AlignCos_FIT_yAllSALL(s,ite,:) = groupmatrixFin.dataFITSALL(s,ite,2) + groupmatrixFin.dataFITSALL(s,ite,3).*(cos(X + groupmatrixFin.dataFITSALL(s,ite,1))); % if you want to plot the predicted function
        end
        subplot(7,6,sp)
        plot(groupmatrixFin.dataALL(s,:)')
        hold on
        plot(groupmatrixFin.AlignCos_FIT_yAllALL(s,:))
        sp=sp+1;
    catch
        probsu(tee) = s;
        tee=tee+1;
        groupmatrixFin.dataFITALL(s,:)=nan;
        groupmatrixFin.dataprefPhaseAll(s)=nan;
    end
    
end

dataAmpSALL      = squeeze(groupmatrixFin.dataFITSALL(:,:,3));
dataAmpOrigALL   = squeeze(groupmatrixFin.dataFITALL(:,3));
dataPhaseOrigALL = squeeze(groupmatrixFin.dataprefPhaseAll);

%compute subject z-score relative to permuted datasets
dataAmpZscoreOrigALL = (dataAmpOrigALL-mean(dataAmpSALL,2))./std(dataAmpSALL,[],2);

save eLifeRev2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%QUESTION 1: FM stimulus modulates behavior in sham?
%basic plottings

Data2Plot = dataGroupPool;
Xnew      = (0:2*pi/1000:2*pi-2*pi/NrAudtACSlag)+2*pi/NrAudtACSlag/2;
Xorig     = 0:2*pi/9:2*pi-2*pi/9;
%colorl    = [246 148 189; 145 130 217; 100 158 166; 209 222 89; 247 195 74; 232 139 93; 27 118 188]./255;
colorl    = [246 148 189; 144 130 217; 100 159 167; 209 222 90; 247 195 75; 232 140 93; 29 110 175]./255;

%%%FIGURE 1: SINGLE SUBJECT EXAMPLES
sp = 1;    
figure

for session = 1:2
    figure
    for subj=1:42%[3 12 42]      
        subplot(6,7,subj)      
        try
        for l =1:6
            plot([Xorig Xorig+2*pi],[Data2Plot(subj,session).hrbyStimbintACSLagMain(l,:) Data2Plot(subj,session).hrbyStimbintACSLagMain(l,:)]','color',colorl(l,:))
            hold on
            fit = Data2Plot(subj,session).hrbyStimbintACSLagMain_FIT(l,:);
            y = fit(2) + fit(3).*(cos(Xnew + fit(1))); % if you want to plot the predicted function
            plot([Xnew Xnew+2*pi],[y y],'--','color',colorl(l,:))
           % scatter(Data2Plot(subj,session).hrbyStimbintACSLagMain_prefPhase(l),max(Data2Plot(subj,session).hrbyStimbintACSLagMain(l,:)))
        end
        plot([Xorig Xorig+2*pi],[Data2Plot(subj,session).hrbyStimbinsham Data2Plot(subj,session).hrbyStimbinsham],'color',colorl(7,:))
        fit = Data2Plot(subj,session).hrbyStimbinsham_FIT;
        y = fit(2) + fit(3).*(cos(Xnew + fit(1))); % if you want to plot the predicted function
        plot([Xnew Xnew+2*pi],[y y],'--','color',colorl(7,:))
        %scatter(Data2Plot(subj,session).hrbyStimbinsham_prefPhase,max(Data2Plot(subj,session).hrbyStimbinsham))
        title(num2str(subj))
        sp = sp+1;
        ylim([0 1.2])
        xtickslabel([])
        catch
        end
    end
end

%%%FIGURE 2: OPTIMAL PHASE
ampGroup           = nan(42,2,7);
PhaseGroup         = nan(42,2,7);
sigGroupSham       = nan(42,2);
sigGroupShamPer95  = nan(length(Data2Plot),2);
sigGroupPerc95AmpS = nan(length(Data2Plot),2);
sigGroupPerc50AmpS = nan(length(Data2Plot),2);
meanAmpSbySession  = nan(length(Data2Plot),2);

for subj=1:length(Data2Plot)
    for session =1:2
        try
            ampGroup(subj,session,:) = [Data2Plot(subj,session).hrbyStimbinsham_FIT(3);Data2Plot(subj,session).hrbyStimbintACSLagMain_FIT(:,3)];
            PhaseGroup(subj,session,:) = [Data2Plot(subj,session).hrbyStimbinsham_prefPhase;Data2Plot(subj,session).hrbyStimbintACSLagMain_prefPhase];
            sigGroupSham(subj,session) = mean(Data2Plot(subj,session).hrbyStimbinsham_FIT(3)<squeeze(Data2Plot(subj,session).hrbyStimbinsham_FITS(:,3)));
            sigGroupShamPer95(subj,session) = prctile(squeeze(Data2Plot(subj,session).hrbyStimbinsham_FITS(:,3)),95);
            sigGroupPerc95AmpS(subj,session) = prctile(squeeze(dataAmpS(subj,:,session)),95);
            sigGroupPerc50AmpS(subj,session) = prctile(squeeze(dataAmpS(subj,:,session)),50);
            meanAmpSbySession(subj,session) = mean(squeeze(dataAmpS(subj,:,session)));            
        catch
        end
    end
end
nv1 = sum(~isnan(sigGroupSham(:,1)));
nv2 = sum(~isnan(sigGroupSham(:,2)));
missingData = ~isnan(ampGroup(:,1,4))&~isnan(ampGroup(:,2,4)); %check lag 4 because there was a participant with not enough trials there

ampGroupsameN           = ampGroup(missingData,:,:);
PhaseGroupsameN         = PhaseGroup(missingData,:,:);
sigGroupShamsameN       = sigGroupSham(missingData,:);
sigGroupShamPerc95sameN = sigGroupShamPer95(missingData,:);

statsAll.nrSinMod2FM_sham = sum(sigGroupShamsameN<0.05);
figure
scatter(ampGroupsameN(:,1,1),ampGroupsameN(:,2,1),'.')
hold on
plot([0 0.5],[0 0.5])
title('Amplitude Sham')
xlabel('S1')
ylabel('S2')
xlim([0 0.5])
ylim([0 0.5])
axis square

for session =1:2
    figure
    for subj=1:length(PhaseGroupsameN)
        subplot(4,10,subj)
        circ_plot(squeeze(PhaseGroupsameN(subj,session,:)),'pretty')
    end
end

figure,
subplot(2,2,1)
distributionPlot(ampGroupsameN(:,:,1), 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
hold on
plot([mean(sigGroupShamPerc95sameN(:,1)) mean(sigGroupShamPerc95sameN(:,1))],[0 3])
plot([mean(sigGroupShamPerc95sameN(:,2)) mean(sigGroupShamPerc95sameN(:,2))],[0 3],'--')

xlim([0 0.5])
title('FM-driven modulation')
ylabel('Session')
xlabel('Amplitude')

subplot(2,2,2)
circ_plot(PhaseGroupsameN(:,1,1),'pretty')
title('S1')
subplot(2,2,3)
circ_plot(PhaseGroupsameN(:,2,1),'pretty')
title('S2')
subplot(2,2,4)
circ_plot(PhaseGroupsameN(:,1,1),'pretty')
title('Inter-session distance')
circ_plot(circ_dist(PhaseGroupsameN(:,1,1),PhaseGroupsameN(:,2,1)),'pretty')

%Q2: is that reliable?
%compare amplitude and phase
[statsAll.S1vsS2ttestamp_H, statsAll.S1vsS2ttestamp_p, statsAll.S1vsS2ttestamp_CI, statsAll.S1vsS2ttestamp_stats]=ttest(ampGroupsameN(:,1,1),ampGroupsameN(:,2,1));
[statsAll.phaseclustShamP(1), statsAll.phaseclustShamZ(1)] = circ_rtest(PhaseGroupsameN(:,1,1));
[statsAll.phaseclustShamP(2), statsAll.phaseclustShamZ(2)] = circ_rtest(PhaseGroupsameN(:,2,1));
[statsAll.phaseclustSessionDistP, statsAll.phaseclustSessionDistV] = circ_vtest(circ_dist(PhaseGroupsameN(:,1,1),PhaseGroupsameN(:,2,1)),0);
[statsAll.phasecorrRho, statsAll.phasecorrP] = circ_corrcc(PhaseGroupsameN(:,1,1),PhaseGroupsameN(:,2,1));
[statsAll.ampcorrRho, statsAll.ampcorrP] = corr(ampGroupsameN(:,1,1),ampGroupsameN(:,2,1),'rows','complete');


%Q3: does tacs modulates optimal stim phase?
figure
for l=1:6
    subplot(2,6,l)
    circ_plot(PhaseGroupsameN(:,1,l+1),'pretty')
    title('S1')
    subplot(2,6,l+6)
    circ_plot(PhaseGroupsameN(:,2,l+1),'pretty')
    title('S1')
end
%check if each phase lag condiftion is clustered around the same mean
for l=1:7
    for session=1:2
        [statsAll.phaseclustlagsP(l,session), statsAll.phaseclustlagsV(l,session)] = circ_vtest(PhaseGroupsameN(:,session,l),circ_mean(PhaseGroupsameN,[],'all'));
    end
end

%Q4: does tACS modulates entrainment amplitude?
%first check if amp across sessions is bigger than pooled
dataAmpOrigbySession   = dataAmpOrig;
dataPhaseOrigbySession = dataPhaseOrig;

figure
subplot(2,3,1)
scatter(mean(dataAmpOrigbySession(missingData,:),2),dataAmpOrigALL(missingData))
ylim([0 0.12])
xlim([0 0.12])
axis square
hold on
plot([0 0.12],[0 0.12])
xlabel('mean across sessions')
ylabel('pooled data')
title('tACS-lag effect amplitude')
subplot(2,3,2)
scatter(dataAmpOrigbySession(missingData,1),dataAmpOrigbySession(missingData,2))
ylim([0 0.12])
xlim([0 0.12])
hold on
plot([0 0.12],[0 0.12])
axis square
xlabel('S1')
ylabel('S2')
title('tACS-lag effect amplitude')
subplot(2,3,3)
distributionPlot(dataAmpOrigbySession(missingData,:), 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
hold on
plot([mean(sigGroupPerc95AmpS(missingData,1)) mean(sigGroupPerc95AmpS(missingData,1)) ],[0 3],'r')
plot([mean(sigGroupPerc95AmpS(missingData,2)) mean(sigGroupPerc95AmpS(missingData,2)) ],[0 3],'r--')
plot([mean(sigGroupPerc50AmpS(missingData,1)) mean(sigGroupPerc50AmpS(missingData,1)) ],[0 3],'k')
plot([mean(sigGroupPerc50AmpS(missingData,2)) mean(sigGroupPerc50AmpS(missingData,2)) ],[0 3],'k--')

xlim([0 0.15])
title('tACS-lag driven modulation')
ylabel('Session')
xlabel('Amplitude')
subplot(2,3,4)
circ_plot(dataPhaseOrigbySession(missingData,1),'pretty')
title('S1')
subplot(2,3,5)
circ_plot(dataPhaseOrigbySession(missingData,2),'pretty')
title('S2')
subplot(2,3,6)
circ_plot(circ_dist(dataPhaseOrigbySession(missingData,1),dataPhaseOrigbySession(missingData,2)),'pretty')
title('Distance')

for session =1:2
    figure
    for subj = 1:42%length(dataAmpOrigbySession(missingData,:))
        subplot(6,7,subj)
        histogram(squeeze(dataAmpS(subj,:,session)))
        hold on
        sigGrouptACSamp(subj,session) = mean(dataAmpS(subj,:,session)>dataAmpOrigbySession(subj,session));
        if mean(dataAmpS(subj,:,session)>dataAmpOrigbySession(subj,session))<0.5
            plot([dataAmpOrigbySession(subj,session) dataAmpOrigbySession(subj,session)],[0 200],'r')
        else
            plot([dataAmpOrigbySession(subj,session) dataAmpOrigbySession(subj,session)],[0 200],'b')
        end
    end
end
statsAll.nrtACSampvsSurr = sum(sigGrouptACSamp(missingData,:)<0.5);

%Q5: are tACS effects reliable across sessions?, optimal tACS lag and
%amplitude
[statsAll.corrtACSAmpS1vsS2_Rho, statsAll.corrtACSAmpS1vsS2_P]=corr(dataAmpOrigbySession(missingData,1),dataAmpOrigbySession(missingData,2),'rows','complete');
[statsAll.corrtACSAmpS1vsS2_ttestH, statsAll.corrtACSAmpS1vsS2_ttestP,statsAll.corrtACSAmpS1vsS2_ttestCI, statsAll.corrtACSAmpS1vsS2_ttestStats]=ttest(dataAmpOrigbySession(missingData,1),dataAmpOrigbySession(missingData,2));
[statsAll.tACSphaseClusterP1, statsAll.tACSphaseClusterZ1]=circ_rtest(dataPhaseOrigbySession(missingData,1));
[statsAll.tACSphaseClusterP2, statsAll.tACSphaseClusterZ2]=circ_rtest(dataPhaseOrigbySession(missingData,2));
[statsAll.tACSphaseClusterdistP, statsAll.tACSphaseClusterdistV]=circ_vtest(circ_dist(dataPhaseOrigbySession(missingData,1),dataPhaseOrigbySession(missingData,2)),0);
[statsAll.tACSphaseCorrRho, statsAll.tACSphaseCorrP]=circ_corrcc(dataPhaseOrigbySession(missingData,1),dataPhaseOrigbySession(missingData,2));
[statsAll.corrtACSAmpperm_ttestH1, statsAll.corrtACSAmpperm_ttestP1,statsAll.corrtACSAmpperm_ttestCI1, statsAll.corrtACSAmpperm_ttestStats1]=ttest(dataAmpOrigbySession(missingData,1),meanAmpSbySession(missingData,1));
[statsAll.corrtACSAmpperm_ttestH2, statsAll.corrtACSAmpperm_ttestP2,statsAll.corrtACSAmpperm_ttestCI2, statsAll.corrtACSAmpperm_ttestStats2]=ttest(dataAmpOrigbySession(missingData,2),meanAmpSbySession(missingData,2));
[statsAll.tACSAmpPooledvsInd_ttestH, statsAll.corrtACSAmpPooledvsInd_ttestP,statsAll.corrtACSAmpPooledvsInd_ttestCI, statsAll.corrtACSAmpPooledvsInd_ttestStats]=ttest(dataAmpOrigALL(missingData),mean(dataAmpOrigbySession(missingData,1),2));

dataAmpOrigbySession2 =dataAmpOrigbySession(missingData,:);
dataPhaseOrigbySession2 =dataPhaseOrigbySession(missingData,:);

% %CHECK difference in gap size?
% [statsAll.corrGapSizeS1vsS2_Rho, statsAll.corrGapSizeS1vsS2_P]=corr(gapSizeGroup(missingData,1),gapSizeGroup(missingData,2),'rows','complete');
% [statsAll.corrGapSizeS1vsS2_ttestH, statsAll.corrGapSizeS1vsS2_ttestP,statsAll.corrGapSizeS1vsS2_ttestCI, statsAll.corrGapSizeS1vsS2_ttestStats]=ttest(gapSizeGroup(missingData,1),gapSizeGroup(missingData,2));
% %correlate change in gap size with change in amp
% [statsAll.corrGapSize2AmpS1vsS2_Rho, statsAll.corrGapSize2AmpS1vsS2_P]=corr(gapSizeGroup(missingData,1)-gapSizeGroup(missingData,2),dataAmpOrigbySession(missingData,1)-dataAmpOrigbySession(missingData,2),'rows','complete');

% figure
% subplot(2,2,1)
% scatter(gapSizeGroup(missingData,1),gapSizeGroup(missingData,2))
% xlim([5 25 ])
% ylim([5 25 ])
% title('GAP size')
% xlabel('Session1')
% ylabel('Session2')
% hold on
% plot([5 25],[5 25])
% axis square
% subplot(2,2,3:4)
% distributionPlot(gapSizeGroup(missingData,:),'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)
% title('GAP size')
% xlabel('Gap size (ms)')
% ylabel('Session')
% hold on
% plot(gapSizeGroup(missingData,:)')
% subplot(2,2,2)
% scatter(zscore(gapSizeGroup(missingData,1)-gapSizeGroup(missingData,2)),zscore(dataAmpOrigbySession(missingData,1)-dataAmpOrigbySession(missingData,2)));
% title('GAP size by tACS amp')
% xlabel('GapSize (zscore)')
% ylabel('tACS amp (zscore)')
% hold on
% plot([-3 3],[-3 3])
% xlim([-3 3])
% ylim([-3 3])
% axis square

%align data using peak or optimal phase from cosine fit
ampReAligSubjCos1   = nan(length(SUBJlist),NrAudtACSlag+1,2);
HRReAligSubjCos1    = nan(length(SUBJlist),NrAudtACSlag+1,2);
HRReAligSubjCos1Max = nan(length(SUBJlist),NrAudtACSlag+1,2);
HRReAligSubjCos1Min = nan(length(SUBJlist),NrAudtACSlag+1,2);
ampReAligSubjRaw    = nan(length(SUBJlist),NrAudtACSlag+1,2);
HRReAligSubjRaw     = nan(length(SUBJlist),NrAudtACSlag+1,2);
HRReAligSubjRawMax  = nan(length(SUBJlist),NrAudtACSlag+1,2);
HRReAligSubjRawMin  = nan(length(SUBJlist),NrAudtACSlag+1,2);
ampReAligSubjCos1S  = nan(length(SUBJlist),NrAudtACSlag+1,NrIte,2);

save eLifeRev2

for session =1:2
    figure
    for subj=1:42
        if ~isnan(dataPhaseOrigbySession(subj,session))
            [optphaselag, optphaselagPOS] = min(abs(circ_dist(lagBins(1:6),dataPhaseOrigbySession(subj,session))));
            [optphaselagRaw, optphaselagPOSRaw] = max(dataGroupPool(subj,session).hrbyStimbintACSLagMain_FIT(:,3));
            ampReAligSubjCos1(subj,:,session) = [dataGroupPool(subj,session).hrbyStimbinsham_FIT(3);wshift('1D',dataGroupPool(subj,session).hrbyStimbintACSLagMain_FIT(:,3),optphaselagPOS-1)];
            HRReAligSubjCos1(subj,:,session) = [dataGroupPool(subj,session).hrbyStimbinsham_FIT(2);wshift('1D',dataGroupPool(subj,session).hrbyStimbintACSLagMain_FIT(:,2),optphaselagPOS-1)];
            HRReAligSubjCos1Max(subj,:,session) = [max(dataGroupPool(subj,session).hrbyStimbinsham); wshift('1D',max(dataGroupPool(subj,session).hrbyStimbintACSLagMain,[],2),optphaselagPOS-1)];
            HRReAligSubjCos1Min(subj,:,session) = [min(dataGroupPool(subj,session).hrbyStimbinsham); wshift('1D',min(dataGroupPool(subj,session).hrbyStimbintACSLagMain,[],2),optphaselagPOS-1)];
            
            ampReAligSubjRaw(subj,:,session) = [dataGroupPool(subj,session).hrbyStimbinsham_FIT(3);wshift('1D',dataGroupPool(subj,session).hrbyStimbintACSLagMain_FIT(:,3),optphaselagPOSRaw-1)];
            HRReAligSubjRaw(subj,:,session) = [dataGroupPool(subj,session).hrbyStimbinsham_FIT(2);wshift('1D',dataGroupPool(subj,session).hrbyStimbintACSLagMain_FIT(:,2),optphaselagPOSRaw-1)];
            HRReAligSubjRawMax(subj,:,session) = [max(dataGroupPool(subj,session).hrbyStimbinsham); wshift('1D',max(dataGroupPool(subj,session).hrbyStimbintACSLagMain,[],2),optphaselagPOSRaw-1)];
            HRReAligSubjRawMin(subj,:,session) = [min(dataGroupPool(subj,session).hrbyStimbinsham); wshift('1D',min(dataGroupPool(subj,session).hrbyStimbintACSLagMain,[],2),optphaselagPOSRaw-1)];
            
            %DO THE SAME FOR THE SURROGATE DATASET
            if sum(dataPhaseOrigS(subj,:,session)) == 0 %subject data is missing
                dataPhaseOrigS(subj,:,session) = nan;
            else
                for ite =1:1000
                    [optphaselagS, optphaselagPOSS] = min(abs(circ_dist(lagBins(1:6),dataPhaseOrigS(subj,ite,session))));
                    ampReAligSubjCos1S(subj,:,ite,session) = [dataGroupPool(subj,session).hrbyStimbinsham_FIT(3);wshift('1D',dataGroupPool(subj,session).hrbyStimbintACSLagMain_FITS(ite,:,3)',optphaselagPOSS-1)];
                end
            end
            
        end
    %plot just to check
    if sum(isnan(ampReAligSubjCos1S(subj,:,ite,session)))> 0
    else
    subplot(6,7,subj)
    plot([dataGroupPool(subj,session).hrbyStimbinsham_FIT(3);dataGroupPool(subj,session).hrbyStimbintACSLagMain_FITS(ite,:,3)'])
    hold on
    plot(ampReAligSubjCos1S(subj,:,ite,session));
    end
    end 
end

figure
subplot(2,4,1)
distributionPlotYCC(squeeze(ampReAligSubjCos1(missingData,:,1)), 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')

subplot(2,4,2)
distributionPlotYCC(squeeze(ampReAligSubjCos1(missingData,:,2)), 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')

subplot(2,4,3)
distributionPlotYCC(squeeze(ampReAligSubjRaw(missingData,:,1)), 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')

subplot(2,4,4)
distributionPlotYCC(squeeze(ampReAligSubjRaw(missingData,:,2)), 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')

subplot(2,4,5)
distributionPlotYCC(squeeze(HRReAligSubjCos1(missingData,:,1)), 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')

subplot(2,4,6)
distributionPlotYCC(squeeze(HRReAligSubjCos1(missingData,:,2)), 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')

subplot(2,4,7)
distributionPlotYCC(squeeze(HRReAligSubjRaw(missingData,:,1)), 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')

subplot(2,4,8)
distributionPlotYCC(squeeze(HRReAligSubjRaw(missingData,:,2)), 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
%correlate change in gap size with change in mod amp (maxadj-minadj)
%subjMod = maxadj-minadj;

maxadj  = squeeze(mean(ampReAligSubjCos1(:,[3 7],:),2));   %average the two bins around the peak
minadj  = squeeze(mean(ampReAligSubjCos1(:,[4 6],:),2));   %average the two bins around the trough
maxadjS = squeeze(mean(ampReAligSubjCos1S(:,[3 7],:,:),2)); %average the two bins around the peak
minadjS = squeeze(mean(ampReAligSubjCos1S(:,[4 6],:,:),2)); %average the two bins around the trough

SubjModNormSham = squeeze(maxadj-minadj)./squeeze(ampReAligSubjCos1(:,1,:));

SubjModAll  = maxadj-minadj;
SubjModAllS = maxadjS-minadjS;

maxvsSham = maxadj-squeeze(ampReAligSubjCos1(:,1,:));
minvsSham = minadj-squeeze(ampReAligSubjCos1(:,1,:));

for ite =1:NrIte
maxvsShamS(:,ite,:) = squeeze(maxadjS(:,ite,:))-squeeze(ampReAligSubjCos1(:,1,:));
minvsShamS(:,ite,:) = squeeze(minadjS(:,ite,:))-squeeze(ampReAligSubjCos1(:,1,:));
end

%transform data to z-score usign the surrogate distributions
subjMAXvsMIN_z  = (SubjModAll-squeeze(mean(SubjModAllS,2)))./squeeze(std(SubjModAllS,[],2));
subjMAXvsSHAM_z = (maxvsSham-squeeze(mean(maxvsShamS,2)))./squeeze(std(maxvsShamS,[],2));
subjMINvsSHAM_z = (minvsSham-squeeze(mean(minvsShamS,2)))./squeeze(std(minvsShamS,[],2));

figure
subplot(1,3,1)
distributionPlotYCC([subjMAXvsSHAM_z(missingData,1) subjMINvsSHAM_z(missingData,1) subjMAXvsMIN_z(missingData,1)], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('Session 1')
ylim([-2.6 2.6])
subplot(1,3,2)
distributionPlotYCC([subjMAXvsSHAM_z(missingData,2) subjMINvsSHAM_z(missingData,2) subjMAXvsMIN_z(missingData,2)], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('Session 2')
ylim([-2.6 2.6])
subplot(1,3,3)
distributionPlotYCC([nanmean(subjMAXvsSHAM_z,2) nanmean(subjMINvsSHAM_z,2) nanmean(subjMAXvsMIN_z,2)], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('average')
ylim([-2.6 2.6])

[statsAll.MAXvsSHAM_z_vs0H1, statsAll.MAXvsSHAM_z_vs0P1, statsAll.MAXvsSHAM_z_vs0CI1, statsAll.MAXvsSHAM_z_vs0Stats1] = ttest(subjMAXvsSHAM_z(:,1));
[statsAll.MINvsSHAM_z_vs0H1, statsAll.MINvsSHAM_z_vs0P1, statsAll.MINXvsSHAM_z_vs0CI1, statsAll.MINvsSHAM_z_vs0Stats1] = ttest(subjMINvsSHAM_z(:,1));
[statsAll.MAXvsMIN_z_vs0H1, statsAll.MAXvsMIN_z_vs0P1, statsAll.MAXvsMIN_z_vs0CI1, statsAll.MAXvsMIN_z_vs0Stats1] = ttest(subjMAXvsMIN_z(:,1));

[statsAll.MAXvsSHAM_z_vs0H2, statsAll.MAXvsSHAM_z_vs0P2, statsAll.MAXvsSHAM_z_vs0CI2, statsAll.MAXvsSHAM_z_vs0Stats2] = ttest(subjMAXvsSHAM_z(:,2));
[statsAll.MINvsSHAM_z_vs0H2, statsAll.MINvsSHAM_z_vs0P2, statsAll.MINXvsSHAM_z_vs0CI2, statsAll.MINvsSHAM_z_vs0Stats2] = ttest(subjMINvsSHAM_z(:,2));
[statsAll.MAXvsMIN_z_vs0H2, statsAll.MAXvsMIN_z_vs0P2, statsAll.MAXvsMIN_z_vs0CI2, statsAll.MAXvsMIN_z_vs0Stats2] = ttest(subjMAXvsMIN_z(:,2));

[statsAll.MAXvsSHAM_z_vs0Hmean, statsAll.MAXvsSHAM_z_vs0Pmean, statsAll.MAXvsSHAM_z_vs0CImean, statsAll.MAXvsSHAM_z_vs0Statsmean] = ttest(nanmean(subjMAXvsSHAM_z,2));
[statsAll.MINvsSHAM_z_vs0Hmean, statsAll.MINvsSHAM_z_vs0Pmean, statsAll.MINXvsSHAM_z_vs0CImean, statsAll.MINvsSHAM_z_vs0Statsmean] = ttest(nanmean(subjMINvsSHAM_z,2));
[statsAll.MAXvsMIN_z_vs0Hmean, statsAll.MAXvsMIN_z_vs0Pmean, statsAll.MAXvsMIN_z_vs0CImean, statsAll.MAXvsMIN_z_vs0Statsmean] = ttest(nanmean(subjMAXvsMIN_z,2));


[r,p]=corr(subjMINvsSHAM_z(missingData,1),subjMINvsSHAM_z(missingData,2))
[r,p]=corr(subjMAXvsSHAM_z(missingData,1),subjMAXvsSHAM_z(missingData,2))
[r,p]=corr(subjMAXvsMIN_z(missingData,1),subjMAXvsMIN_z(missingData,2))



for subj=1:42
    subplot(6,7,subj)
    histogram (maxadjS(subj,:,1))
    hold on
    plot([maxadj(subj,1) maxadj(subj,1)],[0 100])
    
end
for subj=1:42
    subplot(6,7,subj)
    histogram (maxadjS(subj,:,2))
    hold on
    plot([maxadj(subj,2) maxadj(subj,2)],[0 100])
    
end
figure
for subj=1:42
    subplot(6,7,subj)
    histogram (mean(maxadjS(subj,:,2),3))
    hold on
    plot([mean(maxadj(subj,:),2) mean(maxadj(subj,:),2)],[0 100])
    
end

figure
for subj=1:42
    subplot(6,7,subj)
    histogram (mean(SubjModAllS(subj,:,2),3,'omitnan'))
    hold on
    plot([mean(SubjModAll(subj,:),2,'omitnan') mean(SubjModAll(subj,:),2,'omitnan')],[0 100])
    pvalvsSurrog(subj)= mean(mean(SubjModAllS(subj,:,2),3,'omitnan')> mean(SubjModAll(subj,:),2,'omitnan'));
end

%COMPUTE MODULATION Z-SCORE PER SUBJECT USING SURROGATE DATA
ModZscore = (SubjModAll(missingData,:)-squeeze(mean(SubjModAllS(missingData,:,:),2)))./squeeze(std(SubjModAllS(missingData,:,:),[],2));

%DO SOME PLOTTING INCLUDING BOTH REAL AND SURROGATE ANALYSES
figure,
distributionPlotYCC(ModZscore, 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')

figure
subplot(1,2,1)
distributionPlotYCC([ampReAligSubjCos1(missingData,1,1) ampReAligSubjCos1(missingData,7,1) ampReAligSubjCos1(missingData,2:6,1)], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
hold on
plot(mean(squeeze([ampReAligSubjCos1S(missingData,1,1) ampReAligSubjCos1S(missingData,7,1) ampReAligSubjCos1S(missingData,2:6,1)])))
%errorbar(mean(squeeze([ampReAligSubjCos1S(missingData,1,1) ampReAligSubjCos1S(missingData,7,1) ampReAligSubjCos1S(missingData,2:6,1)]))),std(squeeze([ampReAligSubjCos1S(missingData,1,1) ampReAligSubjCos1S(missingData,7,1) ampReAligSubjCos1S(missingData,2:6,1)]))./sqrt(length(ModZscore)));
%
title('FM-driven modulation S1')
ylim([-0.1 0.6])
xlabel('Realigned tACS lags')
xticklabels({'sham' '-pi/3' '0' 'pi/3' '2pi/3' '+-pi' '-2pi/3'})
ylabel('Amplitude')
subplot(1,2,2)
distributionPlotYCC([ampReAligSubjCos1(missingData,1,2) ampReAligSubjCos1(missingData,7,2) ampReAligSubjCos1(missingData,2:6,2)], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('FM-driven modulation S2')
ylim([-0.1 0.6])
xlabel('Realigned tACS lags')
ylabel('Amplitude')
xticklabels({'sham' '-pi/3' '0' 'pi/3' '2pi/3' '+-pi' '-2pi/3'})

figure
subplot(1,2,1)
distributionPlotYCC([squeeze(ampReAligSubjCos1(missingData,1,1)) maxadj(missingData,1) minadj(missingData,1)], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)
ylim([-0.1 0.5])
title('Realigned amplitude S1')
ylabel('Amplitude')
xlabel('Condition')
xticklabels({'sham' 'pos' 'neg'})
subplot(1,2,2)
distributionPlotYCC([squeeze(ampReAligSubjCos1(missingData,1,2)) maxadj(missingData,2) minadj(missingData,2)], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)
ylim([-0.1 0.5])
title('Realigned amplitude S2')
ylabel('Amplitude')
xlabel('Condition')
xticklabels({'sham' 'pos' 'neg'})

figure
distributionPlotYCC(squeeze(ModZscore), 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('tACS-amplitude (Z-score relative to surrogate distribution)')
hold on
subplot(1,2,1)
distributionPlotYCC([squeeze(ampReAligSubjCos1(missingData,1,1)) squeeze(mean(maxadjS(missingData,:,1),2)) squeeze(mean(minadjS(missingData,:,1),2))], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)
ylim([-0.1 0.5])
subplot(1,2,2)
distributionPlotYCC([squeeze(ampReAligSubjCos1(missingData,1,1)) maxadj(missingData,1) minadj(missingData,1)], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)
ylim([-0.1 0.5])
title('Realigned amplitude S1')
ylabel('Amplitude')
xlabel('Condition')
xticklabels({'sham' 'pos' 'neg'})

%[statsAll.corrGapSize2MaxminAdjS1vsS2_Rho, statsAll.corrGapSize2MaxminAdjS1vsS2_P] = corr((gapSizeGroup(missingData,1)-gapSizeGroup(missingData,2)),(SubjModAll(missingData,1)-SubjModAll(missingData,2)),'rows','complete');
[statsAll.corrMaxminAdjS1vsS2_Rho, statsAll.corrMaxminAdjS1vsS2_P] = corr(SubjModAll(missingData,1),SubjModAll(missingData,2),'rows','complete');

figure, scatter(SubjModAll(missingData,1),SubjModAll(missingData,2))
title('Maxadj-Minadj')
xlabel('Session 1')
ylabel('Session 2')

figure, scatter(SubjModAll(missingData,1),SubjModAll(missingData,2))
[rho,p]=corr(SubjModAll(:,1),SubjModAll(:,2),'rows','complete');
[statsAll.ModAdjShamH,statsAll.ModAdjShamP,statsAll.ModAdjShamHCI,statsAll.ModAdjShamStats]=ttest(SubjModAll(missingData,1),SubjModAll(missingData,2));

statsAll.MaxAdjanova_rm_S1=anova_rm([squeeze(ampReAligSubjCos1(missingData,1,1)) maxadj(missingData,1) minadj(missingData,1)]);
statsAll.MaxAdjanova_rm_S2=anova_rm([squeeze(ampReAligSubjCos1(missingData,1,2)) maxadj(missingData,2) minadj(missingData,2)]);
[statsAll.MaxAdjanova_rm_All,statsAll.MaxAdjanova_rm_AllTable]=anova_rm({[squeeze(ampReAligSubjCos1(missingData,1,1)) maxadj(missingData,1) minadj(missingData,1)] [squeeze(ampReAligSubjCos1(missingData,1,2)) maxadj(missingData,2) minadj(missingData,2)]});
[statsAll.posthocttest_S1_H(1), statsAll.posthocttest_S1_P(1), statsAll.posthocttest_S1_CI(1,:), statsAll.posthocttest_S1_Stats(1).stats]=ttest(squeeze(ampReAligSubjCos1(missingData,1,1)),maxadj(missingData,1));
[statsAll.posthocttest_S1_H(2), statsAll.posthocttest_S1_P(2), statsAll.posthocttest_S1_CI(2,:), statsAll.posthocttest_S1_Stats(2).stats]=ttest(squeeze(ampReAligSubjCos1(missingData,1,1)),minadj(missingData,1));
[statsAll.posthocttest_S1_H(3), statsAll.posthocttest_S1_P(3), statsAll.posthocttest_S1_CI(3,:), statsAll.posthocttest_S1_Stats(3).stats]=ttest(maxadj(missingData,1),minadj(missingData,1));
[statsAll.posthocttest_S2_H(1), statsAll.posthocttest_S2_P(1), statsAll.posthocttest_S2_CI(1,:), statsAll.posthocttest_S2_Stats(1).stats]=ttest(squeeze(ampReAligSubjCos1(missingData,1,2)),maxadj(missingData,2));
[statsAll.posthocttest_S2_H(2), statsAll.posthocttest_S2_P(2), statsAll.posthocttest_S2_CI(2,:), statsAll.posthocttest_S2_Stats(2).stats]=ttest(squeeze(ampReAligSubjCos1(missingData,1,2)),minadj(missingData,2));
[statsAll.posthocttest_S2_H(3), statsAll.posthocttest_S2_P(3), statsAll.posthocttest_S2_CI(3,:), statsAll.posthocttest_S2_Stats(3).stats]=ttest(maxadj(missingData,2),minadj(missingData,2));

%%stats and plot for main effect of tacs lag average across sessions
data2test = [mean(ampReAligSubjCos1(missingData,1,1:2),3) mean(maxadj(missingData,:),2) mean(minadj(missingData,:),2)];

[statsAll.posthocttest_ave_H(1), statsAll.posthocttest_ave_P(1), statsAll.posthocttest_ave_CI(1,:), statsAll.posthocttest_ave_Stats(1).stats]=ttest(data2test(:,1),data2test(:,2));
[statsAll.posthocttest_ave_H(2), statsAll.posthocttest_ave_P(2), statsAll.posthocttest_ave_CI(2,:), statsAll.posthocttest_ave_Stats(2).stats]=ttest(data2test(:,1),data2test(:,3));
[statsAll.posthocttest_ave_H(3), statsAll.posthocttest_ave_P(3), statsAll.posthocttest_ave_CI(3,:), statsAll.posthocttest_ave_Stats(3).stats]=ttest(data2test(:,2),data2test(:,3));
statsAll.posthocttest_ave_PBonf = statsAll.posthocttest_ave_P*3;

figure
distributionPlotYCC(data2test, 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
ylim([-0.1 0.5])
title('Realigned amplitude')
ylabel('Amplitude')
xlabel('tACS lag Condition')
xticklabels({'sham' 'pos' 'neg'})

figure
subplot(1,2,1)
distributionPlotYCC([ampReAligSubjCos1(missingData,1,1) ampReAligSubjCos1(missingData,7,1) ampReAligSubjCos1(missingData,2:6,1)], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('FM-driven modulation S1')
ylim([-0.1 0.6])
xlabel('Realigned tACS lags')
xticklabels({'sham' '-pi/3' '0' 'pi/3' '2pi/3' '+-pi' '-2pi/3'})
ylabel('Amplitude')
subplot(1,2,2)
distributionPlotYCC([ampReAligSubjCos1(missingData,1,2) ampReAligSubjCos1(missingData,7,2) ampReAligSubjCos1(missingData,2:6,2)], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('FM-driven modulation S2')
ylim([-0.1 0.6])
xlabel('Realigned tACS lags')
ylabel('Amplitude')
xticklabels({'sham' '-pi/3' '0' 'pi/3' '2pi/3' '+-pi' '-2pi/3'})

figure
subplot(1,2,1)
distributionPlotYCC([squeeze(ampReAligSubjCos1(missingData,1,1)) maxadj(missingData,1) minadj(missingData,1)], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)
ylim([-0.1 0.5])
title('Realigned amplitude S1')
ylabel('Amplitude')
xlabel('Condition')
xticklabels({'sham' 'pos' 'neg'})
subplot(1,2,2)
distributionPlotYCC([squeeze(ampReAligSubjCos1(missingData,1,2)) maxadj(missingData,2) minadj(missingData,2)], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)
ylim([-0.1 0.5])
title('Realigned amplitude S2')
ylabel('Amplitude')
xlabel('Condition')
xticklabels({'sham' 'pos' 'neg'})

[statsAll.InterSessionMod_testH,statsAll.InterSessionMod_testP,~,statsAll.InterSessionMod_testStats]=ttest(SubjModAll(missingData,1),SubjModAll(missingData,2));
[statsAll.Sham2Mod_S1_rho,statsAll.Sham2Mod_S1_p]=corr(ampReAligSubjCos1(missingData,1,1),SubjModAll(missingData,1),'rows','complete');
[statsAll.Sham2Mod_S2_rho,statsAll.Sham2Mod_S2_p]=corr(ampReAligSubjCos1(missingData,1,2),SubjModAll(missingData,2),'rows','complete');

figure
scatter(ampReAligSubjCos1(missingData,1,1),SubjModAll(missingData,1))
hold on
scatter(ampReAligSubjCos1(missingData,1,2),SubjModAll(missingData,2))

% figure, subplot(1,3,1),boxplot(SubjModAll(fami==1,:)),ylim([-1 2.5]),title('fam1')
% subplot(1,3,2),boxplot(SubjModAll(fami==2,:)),ylim([-1 2.5]),title('fam2')
% subplot(1,3,3),boxplot(SubjModAll(fami==3,:)),ylim([-1 2.5]),title('fam3')

%Q6: can we predict tACS effects from Efield, BOLD and Baseline?
%%%%%%%%%%%%%%LOAD EFIELD DATA
load('/mnt/beegfs/users/yuranny.cabral/2019-0226-relentrain/ANA/FEM/Efield_vars_for_GLM_orig_0071_V2_normalPOSNEG.mat')
load('/mnt/beegfs/users/yuranny.cabral/2019-0226-relentrain/ANA/fMRI/GroupData/mriAudiovsBase_Group.mat')

%TAKING MAX-MIN AS MOD
SubjModAll = squeeze(mean(ampReAligSubjCos1(:,[3 7],:),2))-squeeze(mean(ampReAligSubjCos1(:,[4 6],:),2));
%ampCOS2lags2(:,1:2) = dataAmpOrigbySession;%use cos fit on amplit values
ampCOS2lags2(:,1:2) = SubjModAll;% SubjModNormSham;%
ampCOS2lags2(indm,3)=1; %individual montage
ampCOS2lags2(inds,3)=2; %standard montage
ampCOS2lags2(:,4)=1;
ampCOS2lags2(missingData,4)=0; %missing one session
ampCOS2lags2(missingMRI,5)=1; %missing MRI
%newvariable with only participants with MRI
ampCOS2lags3 = ampCOS2lags2((ampCOS2lags2(:,5)==0),:);
ampSham3 = squeeze(ampReAligSubjCos1(ampCOS2lags2(:,5)==0,1,:));
%get Efield variables according to the montage
montage = ampCOS2lags3(:,3);%-Efield2BOLDdist(:,4)
corr2BLD = Efield2BOLDcorr(:,1:2);%;EfieldBOLDdata.corr(:,2:3);
%invert focality and distance
Efield2BOLDdist = -Efield2BOLDdist;
FOC50 = -FOC50;

for mC =1:2
groupEfield(montage==mC,:) = [Efield2BOLDdist(montage==mC,mC), corr2BLD(montage==mC,mC), E_norm_FuncroiMean(montage==mC,mC), E_normal_FuncroiMean(montage==mC,mC),...
    E_norm_RHvsLHratio(montage==mC,mC), E_normal_RHvsLHratio(montage==mC,mC), E_norm(montage==mC,mC), E_normal(montage==mC,mC), FOC50(montage==mC,mC)];
end
ampCOS2lags3 = [ampCOS2lags3 groupEfield];

%remove subjects missing tacs data in MRI
mriAudiovsBase_Group(:,:,:,ampCOS2lags2(:,4)==1)=[];
tempampCos = ampCOS2lags2;
tempampCos(ampCOS2lags2(:,4)==1,:)=[];
tempampCos(tempampCos(:,5)==1,:)=[];

figure,
subplot(1,2,1)
boxplot(ampCOS2lags3(ampCOS2lags3(:,3)==1,1:2))
hold on
plotSpread(ampCOS2lags3(ampCOS2lags3(:,3)==1,1:2))
ylabel('Mod_adj_%Sham')
title('Individual')
xticklabels({'S1' 'S2'})
ylim([-.2 .2])
subplot(1,2,2)
boxplot(ampCOS2lags3(ampCOS2lags3(:,3)==2,1:2))
hold on
plotSpread(ampCOS2lags3(ampCOS2lags3(:,3)==2,1:2))
ylabel('Mod_adj_%Sham')
title('Standard')
ylim([-.2 .2])
xticklabels({'S1' 'S2'})


% %%load beta estimates from fMRI data
load('/mnt/beegfs/users/yuranny.cabral/2019-0226-relentrain/ANA/fMRI/GroupData/groupbetasAudivsBase.mat')
% load('/mnt/beegfs/users/yuranny.cabral/2019-0226-relentrain/ANA/tACS/march2022/PCAresults.mat')
% 
figure
distributionPlot([groupbetas_LH groupbetasRH], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)
title('beta estimates')

%MATRIX ampCos
%Variables: 
%1: tACS effect S1
%2: tACS effect S2 
%3: montage, 
%4: missing one session, 
%5: missing mri, 
%6: Efield2BOLD dist, 
%7: Efield2BOLD corr, 
%8: EnormROI, 
%9: normalROI,
%10: Enorm R-Lratio, 
%11: normal R-L ratio, 
%12: E-norm, 
%13: normal, 
%14: FOC50

%add sham condition to use it as a regressor
ampCOS2lags3 = [ampCOS2lags3 ampSham3 groupbetas_LH groupbetasRH];

%remove subject missing one session
ampCOS2lags3(ampCOS2lags3(:,4)==1,:)=[];
save tACeffectsforSPMana_eLifeRev2 ampCOS2lags3

[statsAll.corrModIndiv_rho,statsAll.corrModIndiv_p]=corr(ampCOS2lags3(ampCOS2lags3(:,3)==1,1),ampCOS2lags3(ampCOS2lags3(:,3)==1,2),'rows','complete');
[statsAll.corrModStand_rho,statsAll.corrModStand_p]=corr(ampCOS2lags3(ampCOS2lags3(:,3)==2,1),ampCOS2lags3(ampCOS2lags3(:,3)==2,2),'rows','complete');

% try partial least square decomposition
predictorsPLS = [zscore(nanmean(ampCOS2lags3(:,1:2),2)) zscore(ampCOS2lags3(:,6)) zscore(ampCOS2lags3(:,7)) zscore(ampCOS2lags3(:,8)) zscore(ampCOS2lags3(:,9)) zscore(ampCOS2lags3(:,12)) zscore(ampCOS2lags3(:,13)) zscore(ampCOS2lags3(:,14)) zscore(nanmean(ampCOS2lags3(:,15:16),2)) zscore(mean(ampCOS2lags3(:,17:18),2))];
X= predictorsPLS(:,2:10);
y= predictorsPLS(:,1);
[XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(X,y,9,'CV',5);

%Plot the percent of variance explained in the response variable (PCTVAR) as a function of the number of components.
figure
plot(1:9,cumsum(100*PCTVAR(2,:)),'-bo');
xlabel('Number of PLS components');
ylabel('Percent Variance Explained in y');

figure
yfit = [ones(size(X,1),1) X]*beta;
residuals = y - yfit;
stem(residuals)
xlabel('Observations');
ylabel('Residuals');

figure, plot(y,yfit,'o')
TSS = sum((y-mean(y)).^2);
RSS = sum((y-yfit).^2);
Rsquared = 1 - RSS/TSS;

%ran the same model and create my surrogate distribution

figure
plot(1:9,stats.W,'o-')
legend({'c1','c2','c3','c4','c5','c6','c7','c8','c9'},'Location','best')
xlabel('Predictor')
ylabel('Weight')

%Calculate the normalized PLS weights.
W0 = stats.W ./ sqrt(sum(stats.W.^2,1));
%Calculate the VIP scores for ncomp components.
p = size(XL,1);
sumSq = sum(XS.^2,1).*sum(yl.^2,1);
vipScore = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));
%Find variables with a VIP score greater than or equal to 1.

indVIP = find(vipScore >= 1);
%Plot the VIP scores.

scatter(1:length(vipScore),vipScore,'x')
hold on
scatter(indVIP,vipScore(indVIP),'rx')
plot([1 length(vipScore)],[1 1],'--k')
hold off
axis tight
xlabel('Predictor Variables')
ylabel('VIP Scores')

collintest(X)

%FIT multiple linear regression models
within = table([1;2],'VariableNames',{'Session'});
between = table(zscore(nanmean(ampCOS2lags3(:,1:2),2)),zscore(ampCOS2lags3(:,6)),zscore(ampCOS2lags3(:,7)),zscore(ampCOS2lags3(:,8)), zscore(ampCOS2lags3(:,9)),...
    zscore(ampCOS2lags3(:,12)),zscore(ampCOS2lags3(:,13)),zscore(ampCOS2lags3(:,14)),...
    zscore(nanmean(ampCOS2lags3(:,15:16),2)),zscore(mean(ampCOS2lags3(:,17:18),2)),...
'VariableNames',{'meanAmp','Dist2Peak','Corr2BOLD','Efield_ROI','Normal_Efield_ROI','Efield','Normal_Efield','Focality','Sham','BOLDbetaMean'});


glm0mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak + Corr2BOLD + Efield_ROI + Normal_Efield_ROI + Efield + Normal_Efield + Focality + Sham + BOLDbetaMean')



glm1mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Efield_ROI + Dist2Peak*Focality + Efield_ROI*Focality')
glm2mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Efield_ROI + Dist2Peak*Focality')
glm3mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Efield_ROI + Efield_ROI*Focality')
glm4mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Focality + Efield_ROI*Focality')
glm5mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Efield_ROI + Dist2Peak*Focality + Efield_ROI*Focality + Sham + BOLDbetaMean')

glm6mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Efield + Dist2Peak*Focality + Efield*Focality')
glm7mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Efield + Dist2Peak*Focality')
glm8mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Efield + Efield*Focality')
glm9mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Focality + Efield*Focality')
glm10mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Efield + Dist2Peak*Focality + Efield*Focality + Sham + BOLDbetaMean')

glm11mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Normal_Efield_ROI + Dist2Peak*Focality + Normal_Efield_ROI*Focality')
glm12mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Normal_Efield_ROI + Dist2Peak*Focality')
glm13mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Normal_Efield_ROI + Normal_Efield_ROI*Focality')
glm14mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Focality + Normal_Efield_ROI*Focality')
glm15mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Normal_Efield_ROI + Dist2Peak*Focality + Normal_Efield_ROI*Focality + Sham + BOLDbetaMean')

glm16mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Normal_Efield + Dist2Peak*Focality + Normal_Efield*Focality')
glm17mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Normal_Efield + Dist2Peak*Focality')
glm18mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Normal_Efield + Normal_Efield*Focality')
glm19mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Focality + Normal_Efield*Focality')
glm20mean = fitlm(between,'meanAmp ~ 1 + Dist2Peak*Normal_Efield + Dist2Peak*Focality + Normal_Efield*Focality + Sham + BOLDbetaMean')

glm21mean = fitlm(between,'meanAmp ~ 1 + Corr2BOLD*Efield_ROI + Efield_ROI*Focality')
glm22mean = fitlm(between,'meanAmp ~ 1 + Corr2BOLD*Efield_ROI + Efield_ROI*Focality + Sham + BOLDbetaMean')

glm23mean = fitlm(between,'meanAmp ~ 1 + Corr2BOLD*Efield + Efield*Focality')
glm24mean = fitlm(between,'meanAmp ~ 1 + Corr2BOLD*Efield')
glm25mean = fitlm(between,'meanAmp ~ 1 + Efield*Focality')
glm26mean = fitlm(between,'meanAmp ~ 1 + Corr2BOLD*Efield + Efield*Focality + Sham + BOLDbetaMean')

glm27mean = fitlm(between,'meanAmp ~ 1 + Corr2BOLD*Normal_Efield_ROI + Normal_Efield_ROI*Focality')
glm28mean = fitlm(between,'meanAmp ~ 1 + Corr2BOLD*Normal_Efield_ROI')
glm29mean = fitlm(between,'meanAmp ~ 1 + Corr2BOLD*Normal_Efield_ROI + Normal_Efield_ROI*Focality + Sham + BOLDbetaMean')

glm30mean = fitlm(between,'meanAmp ~ 1 + Corr2BOLD*Normal_Efield + Normal_Efield*Focality')
glm31mean = fitlm(between,'meanAmp ~ 1 + Corr2BOLD*Normal_Efield')
glm32mean = fitlm(between,'meanAmp ~ 1 + Normal_Efield*Focality')
glm33mean = fitlm(between,'meanAmp ~ 1 + Corr2BOLD*Normal_Efield + Normal_Efield*Focality + Sham + BOLDbetaMean')
glm34mean = fitlm(between,'meanAmp ~ 1 ')

%comparing across montages
AICtableMean = [glm1mean.ModelCriterion.AICc;...
    glm2mean.ModelCriterion.AICc;...
    glm3mean.ModelCriterion.AICc;...
    glm4mean.ModelCriterion.AICc;...
    glm5mean.ModelCriterion.AICc;...
    glm6mean.ModelCriterion.AICc;...
    glm7mean.ModelCriterion.AICc;...
    glm8mean.ModelCriterion.AICc;...
    glm9mean.ModelCriterion.AICc;...
    glm10mean.ModelCriterion.AICc;...
    glm11mean.ModelCriterion.AICc;...
    glm12mean.ModelCriterion.AICc;...
    glm13mean.ModelCriterion.AICc;...
    glm14mean.ModelCriterion.AICc;...
    glm15mean.ModelCriterion.AICc;...
    glm16mean.ModelCriterion.AICc;...
    glm17mean.ModelCriterion.AICc;...
    glm18mean.ModelCriterion.AICc;...
    glm19mean.ModelCriterion.AICc;...
    glm20mean.ModelCriterion.AICc;...
    glm21mean.ModelCriterion.AICc;...
    glm22mean.ModelCriterion.AICc;...
    glm23mean.ModelCriterion.AICc;...
    glm24mean.ModelCriterion.AICc;...
    glm25mean.ModelCriterion.AICc;...
    glm26mean.ModelCriterion.AICc;...
    glm27mean.ModelCriterion.AICc;...
    glm28mean.ModelCriterion.AICc;...
    glm29mean.ModelCriterion.AICc;...
    glm30mean.ModelCriterion.AICc;...
    glm31mean.ModelCriterion.AICc;...
    glm32mean.ModelCriterion.AICc;...
    glm33mean.ModelCriterion.AICc];

[minAIC, minAICpos] = min(AICtableMean);

AICtableMean(:,2)=AICtableMean-AICtableMean(minAICpos);

models = {glm1mean;glm2mean;glm3mean;glm4mean; glm5mean;glm6mean;glm7mean;glm8mean;glm9mean;glm10mean;glm11mean;glm12mean;...
    glm13mean;glm14mean;glm15mean;glm16mean;glm17mean;glm18mean;glm19mean;glm20mean;glm21mean;glm22mean;glm23mean;glm24mean;...
    glm25mean;glm26mean;glm27mean;glm28mean;glm29mean;glm30mean;glm31mean;glm32mean;glm33mean}

%%%%MAKING TABLE
for modelNr =1:length(models)
    ModelTable{modelNr,1} = [ 'tACS_amplitude ~' models{modelNr,1}.Formula.LinearPredictor];
    ModelTable{modelNr,2} = models{modelNr,1}.ModelCriterion.AICc;
    ModelTable{modelNr,3} = models{modelNr,1}.ModelCriterion.AICc-minAIC;
    ModelTable{modelNr,4} = models{modelNr,1}.Rsquared.Ordinary;
    [ModelTable{modelNr,6}, ModelTable{modelNr,5}]= coefTest(models{modelNr,1});
end

ModelTable = sortrows(ModelTable,3);
%check how well the model fits the data
plotSlice(glm18mean)
figure
plotDiagnostics(glm18mean)

% %exclude the two outliers and fit the model again
% glmV2 = fitlm(between,'meanAmp ~ 1 + Enormal*Dist2Peak + Enormal*FOC','Exclude',outlier2Excl');
% plotResiduals(glm19mean2)
% plotSlice(glm19mean2)

% %compare between models with difference in AIC smaller than 2
[h,pValue1,stat1,cValue]= lratiotest(glm21mean.LogLikelihood,glm18mean.LogLikelihood,glm18mean.DFE-glm21mean.DFE);
LR = 2*(glm18mean.LogLikelihood - glm27mean.LogLikelihood); % has a X2 distribution with a df equals to number of constrained parameters, here: 1
pval = 1 - chi2cdf(LR, 1);

LR = 2*(glm18mean.LogLikelihood - glm21mean.LogLikelihood); % has a X2 distribution with a df equals to number of constrained parameters, here: 1
pval = 1 - chi2cdf(LR, 1);

%leave one out cross-validation for the winning model
myArray = table2array(between);
CVO = cvpartition(1:34,'k',34); %number of observations 'n' = 150
err = zeros(CVO.NumTestSets,1);
for i = 1:CVO.NumTestSets
    trIdx = CVO.training(i);
    teIdx = CVO.test(i);
    ytest = classify(meas(teIdx,:),meas(trIdx,:),...
		 species(trIdx,:));
    err(i) = sum(~strcmp(ytest,species(teIdx)));
end
cvErr = sum(err)/sum(CVO.TestSize)



%plots for visualizing the effects
figure
subplot(2,2,1)
scatter(zscore(ampCOS2lags3(:,13)),zscore(mean(ampCOS2lags3(:,1:2),2)))
title('amp by Efield')
xlabel('normal E-field')
ylabel('amp')
xlim([-4 2.5])
ylim([-4 2.5])
axis square

subplot(2,2,2)
scatter(zscore(ampCOS2lags3(:,14)),zscore(mean(ampCOS2lags3(:,1:2),2)))
title('amp by Efield')
xlabel('focality')
ylabel('amp')
xlim([-4 2.5])
ylim([-4 2.5])
axis square
%MATRIX ampCos
%Variables: 
%1: tACS effect S1
%2: tACS effect S2 
%6: Efield2BOLD dist
%13: normal 
%14: FOC50
subplot(2,2,3)
scatter(zscore(ampCOS2lags3(:,6)),zscore(mean(ampCOS2lags3(:,1:2),2)))
title('amp by Efield')
xlabel('dist-to-peak')
ylabel('amp')
xlim([-4 2.5])
ylim([-4 2.5])
axis square 

subplot(2,2,4)
scatter(zscore(ampCOS2lags3(ampCOS2lags3(:,3)==1,13)),zscore(mean(ampCOS2lags3(ampCOS2lags3(:,3)==1,1:2),2)))
hold on
scatter(zscore(ampCOS2lags3(ampCOS2lags3(:,3)==2,13)),zscore(mean(ampCOS2lags3(ampCOS2lags3(:,3)==2,1:2),2)))
title('amp by Efield')
xlabel('normal E-field')
ylabel('amp')
legend({'IO', 'S'})
xlim([-4 2.5])
ylim([-4 2.5])
axis square 

% x = linspace(0,3*pi,200);
% y = cos(x) + rand(1,200);
% sz = 25;
% c = linspace(1,10,length(x));
% scatter(x,y,sz,c,'filled')
% 
% figure
% subplot(1,2,1)
% scatter(normalize(ampCOS2lags3(:,13),'zscore'),normalize(mean(ampCOS2lags3(:,1:2),2),'zscore'),normalize(ampCOS2lags3(:,6),'range')*100+.1,'filled')
% title('Amp by Efield')
% xlabel('normal E-field')
% ylabel('Amp')
% % xlim([-4 2.5])
% % ylim([-4 2.5])
% axis square 
% subplot(1,2,2)
% scatter(normalize(ampCOS2lags3(:,13),'zscore'),normalize(mean(ampCOS2lags3(:,1:2),2),'zscore'),50,normalize(ampCOS2lags3(:,14),'range'),'filled')
% title('Amp by Efield')
% xlabel('normal E-field')
% ylabel('Amp')
% % xlim([-4 2.5])
% % ylim([-4 2.5])
% axis square 
% colorbar

%do median split based on the distance to peak to separate the groups
distGroup = between.Dist2Peak;
distGroup(:,2) = 2;
distGroup(distGroup(:,1)<median(distGroup(:,1)),2) = 1;
distGroup(:,3) = between.Montage;
figure, boxplot([distGroup(distGroup(:,2)==1,1) distGroup(distGroup(:,2)==2,1)])

figure,subplot(2,2,1) %13 14 6
scatter(normalize(ampCOS2lags3(:,13),'zscore'),normalize(mean(ampCOS2lags3(:,1:2),2),'zscore'),normalize(ampCOS2lags3(:,14),'range')*100+.1,normalize(ampCOS2lags3(:,14),'range')*100+.1,'filled')
title('Amp by E-field normal (andFOC)')
xlabel('E-field normal')
ylabel('Amp')
% xlim([-4 2.5])
% ylim([-4 2.5])
hold on
plotInteraction(glm19mean,'FOC','Enormal','predictions') 
axis square 
colorbar
hold on
subplot(2,2,2) %13 14 6
plotInteraction(glm19mean,'FOC','Enormal') 
axis square 
colorbar
hold on

subplot(2,2,3)
scatter(normalize(ampCOS2lags3(:,13),'zscore'),normalize(mean(ampCOS2lags3(:,1:2),2),'zscore'),normalize(ampCOS2lags3(:,6),'range')*100+.1,normalize(ampCOS2lags3(:,6),'range')*100+.1,'filled')
title('Amp by E-field normal (and Dist 2Peak)')
xlabel('E-field normal')
ylabel('Amp')
hold on
hold on
plotInteraction(glm19mean,'Dist2Peak','Enormal','predictions') 
axis square 
colorbar
hold on
subplot(2,2,4)
plotInteraction(glm19mean,'Dist2Peak','Enormal') 
axis square 
colorbar
hold on

%Q7: Intersubject variability
SubjModAllforCV = [squeeze(ampReAligSubjCos1(:,1,:)) maxadj minadj];
SubjModAllforCV(indm,7)=1;
SubjModAllforCV(inds,7)=2;
SubjModAllforCV(missingData,8)=SubjModAllforCV(missingData,7);
SubjModAllforCV(:,9) = SubjModAllforCV(:,3)-SubjModAllforCV(:,5);%submods1
SubjModAllforCV(:,10) = SubjModAllforCV(:,4)-SubjModAllforCV(:,6);%submods2

[statsAll.varh1,statsAll.varp1,statsAll.varci1,statsAll.varstats1] = vartest2(mean(SubjModAllforCV(SubjModAllforCV(:,8)==2,9:10),2),mean(SubjModAllforCV(SubjModAllforCV(:,8)==1,9:10),2));
[statsAll.varhS1,statsAll.varpS1,statsAll.varciS1,statsAll.varstatsS1] = vartest2(SubjModAllforCV(SubjModAllforCV(:,8)==2,9),SubjModAllforCV(SubjModAllforCV(:,8)==1,9));
[statsAll.varhS2,statsAll.varpS2,statsAll.varciS2,statsAll.varstatsS2] = vartest2(SubjModAllforCV(SubjModAllforCV(:,8)==2,10),SubjModAllforCV(SubjModAllforCV(:,8)==1,10));

[statsAll.varhS1S,statsAll.varpS1S,statsAll.varciS1S,statsAll.varstatsS1S] = vartest2(SubjModAllforCV(SubjModAllforCV(:,8)==2,1),SubjModAllforCV(SubjModAllforCV(:,8)==1,1));
[statsAll.varhS2s,statsAll.varpS2s,statsAll.varciS2s,statsAll.varstatsS2s] = vartest2(SubjModAllforCV(SubjModAllforCV(:,8)==2,2),SubjModAllforCV(SubjModAllforCV(:,8)==1,2));
[statsAll.varhS1P,statsAll.varpS1P,statsAll.varciS1P,statsAll.varstatsS1P] = vartest2(SubjModAllforCV(SubjModAllforCV(:,8)==2,3),SubjModAllforCV(SubjModAllforCV(:,8)==1,3));
[statsAll.varhS2P,statsAll.varpS2P,statsAll.varciS2P,statsAll.varstatsS2P] = vartest2(SubjModAllforCV(SubjModAllforCV(:,8)==2,4),SubjModAllforCV(SubjModAllforCV(:,8)==1,4));
[statsAll.varhS1N,statsAll.varpS1N,statsAll.varciS1N,statsAll.varstatsS1N] = vartest2(SubjModAllforCV(SubjModAllforCV(:,8)==2,5),SubjModAllforCV(SubjModAllforCV(:,8)==1,5));
[statsAll.varhS2N,statsAll.varpS2N,statsAll.varciS2N,statsAll.varstatsS2N] = vartest2(SubjModAllforCV(SubjModAllforCV(:,8)==2,6),SubjModAllforCV(SubjModAllforCV(:,8)==1,6));


%compute coefficient of variation for montages
figure
subplot(1,2,1)
distributionPlot(mean(SubjModAllforCV(SubjModAllforCV(:,8)==1,9:10),2), 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)
title('Individual Montage')
ylim([-.3 .4])
subplot(1,2,2)
distributionPlot(mean(SubjModAllforCV(SubjModAllforCV(:,8)==2,9:10),2), 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)
title('Standard Montage')
ylim([-.3 .4])

%COMPUTE COEFFICIENT OF VARIATION
CV(1,1)= std(SubjModAllforCV(SubjModAllforCV(:,8)==1,3))/mean(SubjModAllforCV(SubjModAllforCV(:,8)==1,3));
CV(1,2)= std(SubjModAllforCV(SubjModAllforCV(:,8)==1,4))/mean(SubjModAllforCV(SubjModAllforCV(:,8)==1,4));
CV(1,3)= std(SubjModAllforCV(SubjModAllforCV(:,8)==1,5))/mean(SubjModAllforCV(SubjModAllforCV(:,8)==1,5));
CV(1,4)= std(SubjModAllforCV(SubjModAllforCV(:,8)==1,6))/mean(SubjModAllforCV(SubjModAllforCV(:,8)==1,6));
CV(2,1)= std(SubjModAllforCV(SubjModAllforCV(:,8)==2,3))/mean(SubjModAllforCV(SubjModAllforCV(:,8)==2,3));
CV(2,2)= std(SubjModAllforCV(SubjModAllforCV(:,8)==2,4))/mean(SubjModAllforCV(SubjModAllforCV(:,8)==2,4));
CV(2,3)= std(SubjModAllforCV(SubjModAllforCV(:,8)==2,5))/mean(SubjModAllforCV(SubjModAllforCV(:,8)==2,5));
CV(2,4)= std(SubjModAllforCV(SubjModAllforCV(:,8)==2,6))/mean(SubjModAllforCV(SubjModAllforCV(:,8)==2,6));

%compute iqr

IQR_test(1,1)= iqr(SubjModAllforCV(SubjModAllforCV(:,8)==1,3));
IQR_test(1,2)= iqr(SubjModAllforCV(SubjModAllforCV(:,8)==1,4));
IQR_test(1,3)= iqr(SubjModAllforCV(SubjModAllforCV(:,8)==1,5));
IQR_test(1,4)= iqr(SubjModAllforCV(SubjModAllforCV(:,8)==1,6));
IQR_test(2,1)= iqr(SubjModAllforCV(SubjModAllforCV(:,8)==2,3));
IQR_test(2,2)= iqr(SubjModAllforCV(SubjModAllforCV(:,8)==2,4));
IQR_test(2,3)= iqr(SubjModAllforCV(SubjModAllforCV(:,8)==2,5));
IQR_test(2,4)= iqr(SubjModAllforCV(SubjModAllforCV(:,8)==2,6));


[statsAll.MontageAmph,statsAll.MontageAmpp,statsAll.MontageAmpci,statsAll.MontageAmpstats]=ttest2(ampCOS2lags3(ampCOS2lags3(:,3)==1,1)-ampCOS2lags3(ampCOS2lags3(:,3)==1,2),ampCOS2lags3(ampCOS2lags3(:,3)==2,1)-ampCOS2lags3(ampCOS2lags3(:,3)==2,2));

%investigate predictors for lack of reliability
missingData2 = missingData
missingData2(13) = 0; %exclude outlier in terms of inter session days
ttD = [];
for SUBJ=1:length(SUBJlist)
    if missingData2(SUBJ)==0
    else
        ttD = [ttD; SUBJlist{SUBJ,6:9}];
    end
end

newData = table(zscore(mean(ampCOS2lags2(missingData2,:),2)),zscore(ampCOS2lags2(missingData2,1)-ampCOS2lags2(missingData2,2)), zscore(abs(circ_dist(dataPhaseOrigbySession(missingData2,1),dataPhaseOrigbySession(missingData2,2)))), ttD(:,1), zscore(ttD(:,2)), zscore(ttD(:,3)), zscore(ttD(:,4)), zscore(gapSizeGroup(missingData2,1)-gapSizeGroup(missingData2,2)),ampCOS2lags2(missingData2,3),zscore(mean(gapSizeGroup(missingData2,2),2)),...
'VariableNames',{'meanAmp','distAmp','distCirc','Sex','Age','interDays','interMin','Gapdist','Montage','meanGapSize'});
newData.Sex = categorical(newData.Sex);
newData.Montage = categorical(newData.Montage);

glmDistAmp1 =fitlm(newData,'distAmp ~ 1 + Sex + Age + interDays + interMin + Gapdist + Montage')
glmDistAmp2 =fitlm(newData,'distAmp ~ 1 + Age + interDays + interMin + Gapdist + Montage')
glmDistAmp3 =fitlm(newData,'distAmp ~ 1 + Age + interDays + interMin + Gapdist')
glmDistAmp4 =fitlm(newData,'distAmp ~ 1 + Age + interDays + interMin')
glmDistAmp5 =fitlm(newData,'distAmp ~ 1 + interDays + interMin')
glmDistAmp6 =fitlm(newData,'distAmp ~ 1 + interDays')
glmDistAmp7 =fitlm(newData,'distAmp ~ 1 + interMin')

aicModelDistAmp = [glmDistAmp1.ModelCriterion.AICc; glmDistAmp2.ModelCriterion.AICc;glmDistAmp3.ModelCriterion.AICc;glmDistAmp4.ModelCriterion.AICc;...
    glmDistAmp5.ModelCriterion.AICc;glmDistAmp6.ModelCriterion.AICc; glmDistAmp7.ModelCriterion.AICc];

[h,pValue,stat,cValue] = lratiotest(glmDistAmp4.LogLikelihood,glmDistAmp5.LogLikelihood,1)
[h2,pValue2,stat2,cValue2] = lratiotest(glmDistAmp3.LogLikelihood,glmDistAmp5.LogLikelihood,2)



aicModelDistAmp(:,2)=aicModelDistAmp(:,1)-min(aicModelDistAmp(:,1));
LR = 2*(glmDistAmp4.LogLikelihood - glmDistAmp5.LogLikelihood); % has a X2 distribution with a df equals to number of constrained parameters, here: 1
pval1 = 1 - chi2cdf(LR, 1);
LR = 2*(glmDistAmp3.LogLikelihood - glmDistAmp5.LogLikelihood); % has a X2 distribution with a df equals to number of constrained parameters, here: 1
pval2 = 1 - chi2cdf(LR, 2);

figure,subplot(1,2,1),plotAdjustedResponse(glmDistAmp5,'interDays')
subplot(1,2,2),plotAdjustedResponse(glmDistAmp5,'interMin')

figure,subplot(1,2,1),scatter(ttD(:,3),ampCOS2lags2(missingData2,1)-ampCOS2lags2(missingData2,2))
title('delta Days')
subplot(1,2,2),scatter(ttD(:,4),ampCOS2lags2(missingData2,1)-ampCOS2lags2(missingData2,2))
title('delta Minutes')

%compare between models
glmDistPhase1 =fitlm(newData,'distCirc ~ 1 + Sex + Age + interDays + interMin + Gapdist + Montage')
glmDistPhase2 =fitlm(newData,'distCirc ~ 1 + Sex + Age + interMin + Gapdist + Montage')
glmDistPhase3 =fitlm(newData,'distCirc ~ 1 + Sex + Age + interMin + Gapdist')
glmDistPhase4 =fitlm(newData,'distCirc ~ 1 + Sex + interMin + Gapdist')
glmDistPhase5 =fitlm(newData,'distCirc ~ 1 + Sex + Gapdist')
glmDistPhase6 =fitlm(newData,'distCirc ~ 1 + Sex')

aicModelDistPhase = [glmDistPhase1.ModelCriterion.AICc; glmDistPhase2.ModelCriterion.AICc;glmDistPhase3.ModelCriterion.AICc;glmDistPhase4.ModelCriterion.AICc;...
    glmDistPhase5.ModelCriterion.AICc;glmDistPhase6.ModelCriterion.AICc];

aicModelDistPhase(:,2)=aicModelDistPhase(:,1)-min(aicModelDistPhase(:,1));


save tACSAna_MaxvsMin_0071_v2_eLife_Rev2