
clear
close all
clc

MainFolder = uigetdir([], 'choose main directory with all subjects'); %Main directory with all subjects

subjects = {''};

for subj=1:length(subjects)  
    currSubj = subjects{subj};
    disp(currSubj)
    %cd to subj dir
    cd(fullfile(MainFolder,currSubj))
    
tdcs_lf = sim_struct('TDCSLEADFIELD');
% Head mesh
tdcs_lf.fnamehead = [currSubj '.msh'];
% Output directory
tdcs_lf.pathfem = 'leadfield25';
tdcs_lf.solver_options = 'pardiso';
tdcs_lf.electrode.dimensions = [25,25]; %in mm
tdcs_lf.electrode.thickness  = 2; %in mm
%Define ELECTRODES
%tdcs_lf.electrode = generateElecMontage(66);
run_simnibs(tdcs_lf)
end

%add electrode position between TP7/P7 snd between TP8/P8
for subj=1:length(subjects)
    currSubj = subjects{subj};
    disp(currSubj)
    %cd to subj dir
    cd(fullfile(MainFolder,currSubj))
    %load electrode file
    EEG1010UIJurak2007 = readtable(fullfile(MainFolder,currSubj,['m2m_' currSubj '/eeg_positions/EEG10-10_UI_Jurak_2007.csv']));
    for elec =1:height(EEG1010UIJurak2007)
        if strcmp(EEG1010UIJurak2007.Var5(elec),'TP7')
            posTP7=[EEG1010UIJurak2007.Var2(elec) EEG1010UIJurak2007.Var3(elec) EEG1010UIJurak2007.Var4(elec)];
        elseif strcmp(EEG1010UIJurak2007.Var5(elec),'P7')
            posP7=[EEG1010UIJurak2007.Var2(elec) EEG1010UIJurak2007.Var3(elec) EEG1010UIJurak2007.Var4(elec)];
        elseif strcmp(EEG1010UIJurak2007.Var5(elec),'TP8')
            posTP8=[EEG1010UIJurak2007.Var2(elec) EEG1010UIJurak2007.Var3(elec) EEG1010UIJurak2007.Var4(elec)];
        elseif strcmp(EEG1010UIJurak2007.Var5(elec),'P8')
            posP8=[EEG1010UIJurak2007.Var2(elec) EEG1010UIJurak2007.Var3(elec) EEG1010UIJurak2007.Var4(elec)];
        end
    end
    GroupPos(subj).midTP7_P7 = [mean([posTP7;posP7])];
    GroupPos(subj).midTP8_P8 = [mean([posTP8;posP8])];
end
save([MainFolder '/GroupAna/GroupPos.mat'],'GroupPos')

%---RUN simulation with standard montage
for subj=1:length(subjects)
    currSubj = subjects{subj};
    disp(currSubj)
    %cd to subj dir
    cd(fullfile(MainFolder,currSubj))

S = sim_struct('SESSION'); % Define a stimulation sessions
S.fields = 'veEjJs'; % fields to calculate
S.map_to_fsavg=true;
S.map_to_MNI=true;
S.map_to_surf=true;
S.map_to_vol=true;
S.fnamehead = [currSubj '.msh']; % Choose the head mesh

%run the simulation for several positions
S.pathfem = 'simulations_FC6_TP8-P8_FC5_TP7-P7/'; % Folder for the simulation output

%Montage from Baltus...Herrmann 2018 FC5-TP7/P7 Chan3: FC6-TP8/P8
% Define the tDCS simulation
S.poslist{1} = sim_struct('TDCSLIST');
S.poslist{1}.currents = [0.5e-3, -0.5e-3, 0.5e-3, -0.5e-3]; % Current going through each channel, in Ampere

%First Electrode (RH1)
S.poslist{1}.electrode(1).channelnr = 1; % Connect it to the fisrt channel (1 mA)
S.poslist{1}.electrode(1).centre = 'FC6';
S.poslist{1}.electrode(1).shape = 'ellipse'; % Define it to be elliptical/circular
S.poslist{1}.electrode(1).dimensions = [25, 25]; % Electrode diameter (in mm)
S.poslist{1}.electrode(1).thickness = 2; % Electrode 

%Second Electrode (RH2)
S.poslist{1}.electrode(2).channelnr = 2; %Connect it to the second channel (-1mA)
%S.poslist{1}.electrode(2).centre = 'TP8'; % Place it over T7
S.poslist{1}.electrode(2).centre = GroupPos(subj).midTP8_P8;
S.poslist{1}.electrode(2).shape = 'ellipse'; % Elliptical shape
S.poslist{1}.electrode(2).dimensions = [25, 25]; %Diameter of 25mm
S.poslist{1}.electrode(2).thickness = 2;

%First Electrode (LH1)
S.poslist{1}.electrode(3).channelnr = 3; % Connect it to the fisrt channel (1 mA)
S.poslist{1}.electrode(3).centre = 'FC5'; % Place it over T8
S.poslist{1}.electrode(3).shape = 'ellipse'; % Define it to be elliptical/circular
S.poslist{1}.electrode(3).dimensions = [25, 25]; % Electrode diameter (in mm)
S.poslist{1}.electrode(3).thickness = 2; % Electrode 

%Second Electrode (LH2)
S.poslist{1}.electrode(4).channelnr = 4; %Connect it to the second channel (-1mA)
%S.poslist{1}.electrode(4).centre = 'TP7'; % Place it over T8
S.poslist{1}.electrode(4).centre = GroupPos(subj).midTP7_P7;
S.poslist{1}.electrode(4).shape = 'ellipse'; % Elliptical shape
S.poslist{1}.electrode(4).dimensions = [25, 25]; %Diameter of 25mm
S.poslist{1}.electrode(4).thickness = 2;

% Run Simulation
run_simnibs(S);
end
%%run zoefel's 2021 PlosBiol montage
for subj=1:length(subjects)
    currSubj = subjects{subj};
    disp(currSubj)
    %cd to subj dir
    cd(fullfile(MainFolder,currSubj))

S = sim_struct('SESSION'); % Define a stimulation sessions
S.fields = 'veEjJs'; % fields to calculate
S.map_to_fsavg=true;
S.map_to_MNI=true;
S.map_to_surf=true;
S.map_to_vol=true;
S.fnamehead = [currSubj '.msh']; % Choose the head mesh

%run the simulation for several positions
S.pathfem = 'simulations_T7-T8_ring2/'; % Folder for the simulation output

%Montage from Baltus...Herrmann 2018 FC5-TP7/P7 Chan3: FC6-TP8/P8
% Define the tDCS simulation
S.poslist{1} = sim_struct('TDCSLIST');
%S.poslist{1}.currents = [0.5e-3, -0.5e-3, 0.5e-3, -0.5e-3]; % Current going through each channel, in Ampere
S.poslist{1}.currents = [-0.5e-3, 0.5e-3, -0.5e-3, 0.5e-3]; % Current going through each channel, in Ampere

%% Define the tDCS simulation
S.poslist{1} = sim_struct('TDCSLIST');
%S.poslist{1}.currents = [0.5e-3, -0.5e-3, 0.5e-3, -0.5e-3]; % Current going through each channel, in Ampere
S.poslist{1}.currents = [-0.5e-3, 0.5e-3, -0.5e-3, 0.5e-3]; % Current going through each channel, in Ampere

%First Electrode (circular, central)
S.poslist{1}.electrode(1).channelnr = 1; % Connect it to the fisrt channel (1 mA)
S.poslist{1}.electrode(1).centre = 'T8'; % Place it over T7
S.poslist{1}.electrode(1).shape = 'ellipse'; % Define it to be elliptical/circular
S.poslist{1}.electrode(1).dimensions = [20, 20]; % Electrode diameter (in mm)
S.poslist{1}.electrode(1).thickness = 2; % Electrode 

%Second Electrode (ring, surrounding)
S.poslist{1}.electrode(2).channelnr = 2; %Connect it to the second channel (-1mA)
S.poslist{1}.electrode(2).centre = 'T8'; % Place it over T7
S.poslist{1}.electrode(2).shape = 'ellipse'; % Elliptical shape
S.poslist{1}.electrode(2).dimensions = [100, 100]; %Diameter of 100mm
S.poslist{1}.electrode(2).thickness = 2;

% firstHole
S.poslist{1}.electrode(2).holes = sim_struct('ELECTRODE'); % Define the hole
S.poslist{1}.electrode(2).holes.centre = 'T8'; % Hole is also centered in T7
S.poslist{1}.electrode(2).holes.shape = 'ellipse'; % Shape of the hole
S.poslist{1}.electrode(2).holes.dimensions = [75, 75]; % Diameter of 75mm

%First Electrode (circular, central)
S.poslist{1}.electrode(3).channelnr = 3; % Connect it to the fisrt channel (1 mA)
S.poslist{1}.electrode(3).centre = 'T7'; % Place it over T8
S.poslist{1}.electrode(3).shape = 'ellipse'; % Define it to be elliptical/circular
S.poslist{1}.electrode(3).dimensions = [20, 20]; % Electrode diameter (in mm)
S.poslist{1}.electrode(3).thickness = 2; % Electrode 

%Second Electrode (ring, surrounding)
S.poslist{1}.electrode(4).channelnr = 4; %Connect it to the second channel (-1mA)
S.poslist{1}.electrode(4).centre = 'T7'; % Place it over T8
S.poslist{1}.electrode(4).shape = 'ellipse'; % Elliptical shape
S.poslist{1}.electrode(4).dimensions = [100, 100]; %Diameter of 100mm
S.poslist{1}.electrode(4).thickness = 2;

% second Hole
S.poslist{1}.electrode(4).holes = sim_struct('ELECTRODE'); % Define the hole
S.poslist{1}.electrode(4).holes.centre = 'T7'; % Hole is also centered in T7
S.poslist{1}.electrode(4).holes.shape = 'ellipse'; % Shape of the hole
S.poslist{1}.electrode(4).holes.dimensions = [75, 75]; % Diameter of 75mm


% Run Simulation
run_simnibs(S);
end
%run max stim to get max achievable target intensity
%load ROI centroid coordinate
%targetCoord = getIndivTargetCoordinates;
load(fullfile('ROIcentroid')
for subj=1:length(subjects)
    currSubj = subjects{subj};
    disp(currSubj)
    if ~strcmp(currSubj,ROIcentroid(subj).ID) %check that the subject IDs match
        error 'SUbject ID do not match'
    else
        %get left and right ROI coordinates
        if ROIcentroid(subj).roi(1,1)<0 %first roi is left hemisphere
           targetCoord(subj).LH = ROIcentroid(subj).roi(1,:);
           targetCoord(subj).RH = ROIcentroid(subj).roi(2,:);
        else
           targetCoord(subj).LH = ROIcentroid(subj).roi(2,:);
           targetCoord(subj).RH = ROIcentroid(subj).roi(1,:);
        end
    %cd to subj dir
    cd(fullfile(MainFolder,currSubj))
%read the mesh
mesh_load_hdf5(fullfile('leadfield25/',[currSubj '_leadfield_EEG10-10_UI_Jurak_2007.hdf5']))

%2.-OPTIMIZATION functional ROI
% optimization with two targets
opt               = opt_struct('TDCSoptimize');
opt.leadfield_hdf = fullfile('leadfield25/',[currSubj '_leadfield_EEG10-10_UI_Jurak_2007.hdf5']);
opt.name          = 'optimization_Func_Step1_Orient0071/two_targetsFinalMax';

opt.max_total_current      = 2e-3;
opt.max_individual_current = 0.5e-3;
opt.max_active_electrodes  = 4;
% Target in the left auditory cortex
opt.target(1).positions  = mni2subject_coords(targetCoord(subj).LH, ['m2m_' currSubj]);
%opt.target(1).positions  = mni2subject_coords([-57 -22 5], ['m2m_' currSubj]);
opt.target(1).intensity  = 0.1; %to maximize intensity in the target
opt.target(1).type       = 'TDCStarget';
opt.target(1).directions = [0, 0.7, 1];%'normal';%
opt.target(1).radius     = 3;

% Target in the right auditory cortex
opt.target(2).positions  = mni2subject_coords(targetCoord(subj).RH, ['m2m_' currSubj]);
%opt.target(2).positions  = mni2subject_coords([57 -13 5], ['m2m_' currSubj]);
opt.target(2).intensity  = 0.1; %to maximize intensity in the target
opt.target(2).type       = 'TDCStarget';
opt.target(2).directions = [0, 0.7, 1];%'normal';
opt.target(2).radius     = 3;
opt.avoid(1).tissues     = 1006; % 1006 corresponds to the eye surface
% opt.avoid(2).positions = opt.target(2).positions; % position of target receiving more stimulation
% opt.avoid(2).weights = 1; % Do not set this too large
% opt.avoid(3).positions = opt.target(2).positions; % position of target receiving more stimulation
% opt.avoid(3).weights = 1; % Do not set this too large
run_simnibs(opt)
    end
end
save('targetCoord','targetCoord')

% targetIntensity = getTargetIntensityFunct2; %load target intensity
%save('targetintensity2submit','targetIntensity')

load('Intensity_optimization_Func_Step1_Orient0071_MAXminus01','targetIntensity')
load('targetCoord','targetCoord')%load traget coordinates

for subj= 1:length(subjects)
    currSubj = subjects{subj};
    disp(currSubj)
     if ~strcmp(currSubj,targetIntensity(subj).ID) %check that the subject IDs match
        error 'SUbject ID do not match'
    else
    %cd to subj dir
    cd(fullfile(MainFolder,currSubj))
%read the mesh
mesh_load_hdf5(fullfile('leadfield25/',[currSubj '_leadfield_EEG10-10_UI_Jurak_2007.hdf5']))

%2.-OPTIMIZATION
% optimization with two targets
opt               = opt_struct('TDCSoptimize');
opt.leadfield_hdf = fullfile('leadfield25/',[currSubj '_leadfield_EEG10-10_UI_Jurak_2007.hdf5']);
opt.name          = 'optimization_Func_Step1_Orient0071_maxminus01/two_targetsFinal';

opt.max_total_current      = 2e-3;
opt.max_individual_current = 0.5e-3;
opt.max_active_electrodes  = 4;

% Target in the right auditory cortex
opt.target(1).positions  = mni2subject_coords(targetCoord(subj).LH, ['m2m_' currSubj]);
%opt.target(1).positions  = mni2subject_coords([-57 -22 5], ['m2m_' currSubj]);
opt.target(1).intensity  = targetIntensity(subj).Int; %to maximize intensity in the target
opt.target(1).type       = 'TDCStarget';
opt.target(1).directions = [0, 0.7, 1];%'normal';%
opt.target(1).radius     = 3;

% Target in the left auditory cortex
opt.target(2).positions  = mni2subject_coords(targetCoord(subj).RH, ['m2m_' currSubj]);
%opt.target(2).positions  = mni2subject_coords([57 -13 5], ['m2m_' currSubj]);
opt.target(2).intensity  = targetIntensity(subj).Int; %to maximize intensity in the target
opt.target(2).type       = 'TDCStarget';
opt.target(2).directions = [0, 0.7, 1];%'normal';
opt.target(2).radius     = 3;
opt.avoid(1).tissues     = 1006; % 1006 corresponds to the eye surface
% opt.avoid(2).positions = opt.target(2).positions; % position of target receiving more stimulation
% opt.avoid(2).weights = 1; % Do not set this too large
% opt.avoid(3).positions = opt.target(2).positions; % position of target receiving more stimulation
% opt.avoid(3).weights = 1; % Do not set this too large

run_simnibs(opt)
     end
end
%REPEAT OPTIMIZATION FOR PARTICIPANTS THAT DID NOT REACH TARGET STRENGTH IN
%BOTH TARGETS OR RESULTS SHOW TOO DIFFFERENT FIELDS BETWEEN ELECTRODES
subjects2 = [3 6 20 25 33 38];
load(fullfile('Intensity_optimization_Func_Step1_Orient0071_MAXminus01','targetIntensity')
load(fullfile('targetCoord','targetCoord')%load traget coordinates

for subj= subjects2
    currSubj = subjects{subj};
    disp(currSubj)
     if ~strcmp(currSubj,targetIntensity(subj).ID) %check that the subject IDs match
        error 'SUbject ID do not match'
     else
        targetIntensity(subj).Int = targetIntensity(subj).Int-0.01;
    %cd to subj dir
    cd(fullfile(MainFolder,currSubj))
%read the mesh
mesh_load_hdf5(fullfile('leadfield25/',[currSubj '_leadfield_EEG10-10_UI_Jurak_2007.hdf5']))

%2.-OPTIMIZATION
% optimization with two targets
opt               = opt_struct('TDCSoptimize');
opt.leadfield_hdf = fullfile('leadfield25/',[currSubj '_leadfield_EEG10-10_UI_Jurak_2007.hdf5']);
opt.name          = 'optimization_Func_Step2_Orient0071_maxminus02/two_targetsFinal';

opt.max_total_current      = 2e-3;
opt.max_individual_current = 0.5e-3;
opt.max_active_electrodes  = 4;

% Target in the right auditory cortex
opt.target(1).positions  = mni2subject_coords(targetCoord(subj).LH, ['m2m_' currSubj]);
%opt.target(1).positions  = mni2subject_coords([-57 -22 5], ['m2m_' currSubj]);
opt.target(1).intensity  = targetIntensity(subj).Int; %to maximize intensity in the target
opt.target(1).type       = 'TDCStarget';
opt.target(1).directions = [0, 0.7, 1];%'normal';%
opt.target(1).radius     = 3;

% Target in the left auditory cortex
opt.target(2).positions  = mni2subject_coords(targetCoord(subj).RH, ['m2m_' currSubj]);
%opt.target(2).positions  = mni2subject_coords([57 -13 5], ['m2m_' currSubj]);
opt.target(2).intensity  = targetIntensity(subj).Int; %to maximize intensity in the target
opt.target(2).type       = 'TDCStarget';
opt.target(2).directions = [0, 0.7, 1];%'normal';
opt.target(2).radius     = 3;
opt.avoid(1).tissues     = 1006; % 1006 corresponds to the eye surface
% opt.avoid(2).positions = opt.target(2).positions; % position of target receiving more stimulation
% opt.avoid(2).weights = 1; % Do not set this too large
% opt.avoid(3).positions = opt.target(2).positions; % position of target receiving more stimulation
% opt.avoid(3).weights = 1; % Do not set this too large

run_simnibs(opt)
     end
end

%REPEAT OPTIMIZATION FOR PARTICIPANTS THAT DID NOT REACH TARGET STRENGTH IN
%BOTH TARGETS OR RESULTS SHOW TOO DIFFFERENT FIELDS BETWEEN ELECTRODES
subjects2 = [25 33];
load(fullfile('Intensity_optimization_Func_Step1_Orient0071_MAXminus01','targetIntensity')
load(fullfile('targetCoord','targetCoord')%load target coordinates

for subj= subjects2
    currSubj = subjects{subj};
    disp(currSubj)
     if ~strcmp(currSubj,targetIntensity(subj).ID) %check that the subject IDs match
        error 'SUbject ID do not match'
     else
        targetIntensity(subj).Int = targetIntensity(subj).Int-0.02;
    %cd to subj dir
    cd(fullfile(MainFolder,currSubj))
%read the mesh
mesh_load_hdf5(fullfile('leadfield25/',[currSubj '_leadfield_EEG10-10_UI_Jurak_2007.hdf5']))

%2.-OPTIMIZATION
% optimization with two targets
opt               = opt_struct('TDCSoptimize');
opt.leadfield_hdf = fullfile('leadfield25/',[currSubj '_leadfield_EEG10-10_UI_Jurak_2007.hdf5']);
opt.name          = 'optimization_Func_Step3_Orient0071_maxminus03/two_targetsFinal';

opt.max_total_current      = 2e-3;
opt.max_individual_current = 0.5e-3;
opt.max_active_electrodes  = 4;

% Target in the right auditory cortex
opt.target(1).positions  = mni2subject_coords(targetCoord(subj).LH, ['m2m_' currSubj]);
%opt.target(1).positions  = mni2subject_coords([-57 -22 5], ['m2m_' currSubj]);
opt.target(1).intensity  = targetIntensity(subj).Int; %to maximize intensity in the target
opt.target(1).type       = 'TDCStarget';
opt.target(1).directions = [0, 0.7, 1];%'normal';%
opt.target(1).radius     = 3;

% Target in the left auditory cortex
opt.target(2).positions  = mni2subject_coords(targetCoord(subj).RH, ['m2m_' currSubj]);
%opt.target(2).positions  = mni2subject_coords([57 -13 5], ['m2m_' currSubj]);
opt.target(2).intensity  = targetIntensity(subj).Int; %to maximize intensity in the target
opt.target(2).type       = 'TDCStarget';
opt.target(2).directions = [0, 0.7, 1];%'normal';
opt.target(2).radius     = 3;
opt.avoid(1).tissues     = 1006; % 1006 corresponds to the eye surface
% opt.avoid(2).positions = opt.target(2).positions; % position of target receiving more stimulation
% opt.avoid(2).weights = 1; % Do not set this too large
% opt.avoid(3).positions = opt.target(2).positions; % position of target receiving more stimulation
% opt.avoid(3).weights = 1; % Do not set this too large

run_simnibs(opt)
     end
end
        
%Efield simulations with winning montage
%opMontage = getOptMontage2;
opMontage = getOptMontage_0071;
for subj=1:length(subjects)
    currSubj = subjects{subj};
    disp(currSubj)
    %cd to subj dir
    cd(fullfile(MainFolder,currSubj))

S = sim_struct('SESSION'); % Define a stimulation sessions
S.fields = 'veEjJs'; % fields to calculate
S.map_to_fsavg=true;
S.map_to_MNI=true;
S.map_to_surf=true;
S.map_to_vol=true;
S.fnamehead = [currSubj '.msh']; % Choose the head mesh

% run the simulation for several positions
S.pathfem = 'simulations_Func_0071_V2/'; % Folder for the simulation output


%% Define the tDCS simulation
S.poslist{1} = sim_struct('TDCSLIST');
S.poslist{1}.currents = [0.5e-3, -0.5e-3, 0.5e-3, -0.5e-3]; % Current going through each channel, in Ampere

%First Electrode (RH1)
S.poslist{1}.electrode(1).channelnr = 1; % Connect it to the fisrt channel (1 mA)
S.poslist{1}.electrode(1).centre = opMontage(subj).RH{1};
S.poslist{1}.electrode(1).shape = 'ellipse'; % Define it to be elliptical/circular
S.poslist{1}.electrode(1).dimensions = [25, 25]; % Electrode diameter (in mm)
S.poslist{1}.electrode(1).thickness = 2; % Electrode 

%Second Electrode (RH2)
S.poslist{1}.electrode(2).channelnr = 2; %Connect it to the second channel (-1mA)
S.poslist{1}.electrode(2).centre = opMontage(subj).RH{2}; % Place it over T7
S.poslist{1}.electrode(2).shape = 'ellipse'; % Elliptical shape
S.poslist{1}.electrode(2).dimensions = [25, 25]; %Diameter of 25mm
S.poslist{1}.electrode(2).thickness = 2;

%First Electrode (LH1)
S.poslist{1}.electrode(3).channelnr = 3; % Connect it to the first channel (1 mA)
S.poslist{1}.electrode(3).centre = opMontage(subj).LH{1}; % Place it over T8
S.poslist{1}.electrode(3).shape = 'ellipse'; % Define it to be elliptical/circular
S.poslist{1}.electrode(3).dimensions = [25, 25]; % Electrode diameter (in mm)
S.poslist{1}.electrode(3).thickness = 2; % Electrode 

%Second Electrode (LH2)
S.poslist{1}.electrode(4).channelnr = 4; %Connect it to the second channel (-1mA)
S.poslist{1}.electrode(4).centre = opMontage(subj).LH{2}; % Place it over T8
S.poslist{1}.electrode(4).shape = 'ellipse'; % Elliptical shape
S.poslist{1}.electrode(4).dimensions = [25, 25]; %Diameter of 25mm
S.poslist{1}.electrode(4).thickness = 2;

%% Run Simulation
run_simnibs(S);
end
%

