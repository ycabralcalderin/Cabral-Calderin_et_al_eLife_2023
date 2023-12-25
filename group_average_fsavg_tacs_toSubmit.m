%% Load simulation results
subjects = {''};
results_folder = fullfile('simulations_FC6_TP8-P8_FC5_TP7-P7', 'fsavg_overlays');
fsavg_msh_name = '_TDCS_1_scalar_fsavg.msh';
field_name = 'E_normal';
fields = {};
for i = 1:length(subjects)
    sub = subjects{i};
    % load mesh with results transformed to fsaverage space
    m = mesh_load_gmsh4(fullfile(pwd, sub, results_folder, [sub fsavg_msh_name]));
    % Save the field of each subject
    fields{i} = m.node_data{get_field_idx(m, field_name, 'node')}.data;
end
%% Calculate and plot averages
% Calculate
fields = cell2mat(fields);
avg_field = mean(fields, 2);
std_field = std(fields, 0, 2);
% Plot
m.node_data = {}; %cleanup fields
m.node_data{1}.data = avg_field; % add average field
m.node_data{1}.name = [field_name '_avg'];
m.node_data{2}.data = std_field; % add std field
m.node_data{2}.name = [field_name '_std'];

% show surfaces with fields
h1 = mesh_show_surfaceycc(m, 'field_idx', [field_name '_avg'],'scaleLimits',[-0.1 0.1]);
title('E_normal_avg Standard')
h2 = mesh_show_surfaceycc(m, 'field_idx', [field_name '_std'],'scaleLimits',[-0.1 0.1]);
title('E_normal_std Standard')
savefig(h1,'E_normal_avg Standard')
savefig(h2,'E_normal_std Standard')
save('E_normal_avg_Standard','m')
%%%%%%%%%FILED ENORM
field_name = 'E_norm';
fields = {};
for i = 1:length(subjects)
    sub = subjects{i};
    % load mesh with results transformed to fsaverage space
    m = mesh_load_gmsh4(fullfile(pwd, sub, results_folder, [sub fsavg_msh_name]));
    % Save the field of each subject
    fields{i} = m.node_data{get_field_idx(m, field_name, 'node')}.data;
end
%% Calculate and plot averages
% Calculate
fields = cell2mat(fields);
avg_field = mean(fields, 2);
std_field = std(fields, 0, 2);
% Plot
m.node_data = {}; %cleanup fields
m.node_data{1}.data = avg_field; % add average field
m.node_data{1}.name = [field_name '_avg'];
m.node_data{2}.data = std_field; % add std field
m.node_data{2}.name = [field_name '_std'];

% show surfaces with fields
h3=mesh_show_surfaceycc(m, 'field_idx', [field_name '_avg'],'scaleLimits',[-0.19 0.19]);
title('E_norm_avg Standard')
h4=mesh_show_surfaceycc(m, 'field_idx', [field_name '_std'],'scaleLimits',[-0.1 0.1]);
title('E_norm_std Standard')
savefig(h3,'E_norm_avg Standard')
savefig(h4,'E_norm_std Standard')
save('E_norm_avg_Standard','m')

%%%%%____________
%%OPTIMIZED MONTAGE
results_folder = fullfile('simulations_Func_0071_V2', 'fsavg_overlays');
fsavg_msh_name = '_TDCS_1_scalar_fsavg.msh';
field_name = 'E_normal';
fields = {};
for i = 1:length(subjects)
    sub = subjects{i};
    % load mesh with results transformed to fsaverage space
    m = mesh_load_gmsh4(fullfile(pwd, sub, results_folder, [sub fsavg_msh_name]));
    % Save the field of each subject
    fields{i} = m.node_data{get_field_idx(m, field_name, 'node')}.data;
end
%% Calculate and plot averages
% Calculate
fields = cell2mat(fields);
avg_field = mean(fields, 2);
std_field = std(fields, 0, 2);
% Plot
m.node_data = {}; %cleanup fields
m.node_data{1}.data = avg_field; % add average field
m.node_data{1}.name = [field_name '_avg'];
m.node_data{2}.data = std_field; % add std field
m.node_data{2}.name = [field_name '_std'];

% show surfaces with fields
h5=mesh_show_surfaceycc(m, 'field_idx', [field_name '_avg'],'scaleLimits',[-0.1 0.1]);
title('E_normal_avg Optimized')
h6=mesh_show_surfaceycc(m, 'field_idx', [field_name '_std'],'scaleLimits',[-0.1 0.1]);
title('E_normal_std Optimized')
savefig(h5,'E_normal_avg Optimized')
savefig(h6,'E_normal_std Optimized')
save('E_normal_avg_Optimized','m')

%%%%%%%%%FILED ENORM
field_name = 'E_norm';
fields = {};
for i = 1:length(subjects)
    sub = subjects{i};
    % load mesh with results transformed to fsaverage space
    m = mesh_load_gmsh4(fullfile(pwd, sub, results_folder, [sub fsavg_msh_name]));
    % Save the field of each subject
    fields{i} = m.node_data{get_field_idx(m, field_name, 'node')}.data;
end
%% Calculate and plot averages
% Calculate
fields = cell2mat(fields);
avg_field = mean(fields, 2);
std_field = std(fields, 0, 2);
% Plot
m.node_data = {}; %cleanup fields
m.node_data{1}.data = avg_field; % add average field
m.node_data{1}.name = [field_name '_avg'];
m.node_data{2}.data = std_field; % add std field
m.node_data{2}.name = [field_name '_std'];

% show surfaces with fields
h7=mesh_show_surfaceycc(m, 'field_idx', [field_name '_avg'],'scaleLimits',[-0.19 0.19]);
title('E_norm_avg Optimized')
h8=mesh_show_surfaceycc(m, 'field_idx', [field_name '_std'],'scaleLimits',[-0.1 0.1]);
title('E_norm_std Optimized')
savefig(h7,'E_norm_avg Optimized')
savefig(h8,'E_norm_std Optimized')
save('E_norm_avg_Optimized','m')

%Elec Ring
results_folder = fullfile('simulations_T7-T8_ring2', 'fsavg_overlays');
fsavg_msh_name = '_TDCS_1_scalar_fsavg.msh';
field_name = 'E_normal';
fields = {};
for i = 1:length(subjects)
    sub = subjects{i};
    % load mesh with results transformed to fsaverage space
    m = mesh_load_gmsh4(fullfile(pwd, sub, results_folder, [sub fsavg_msh_name]));
    % Save the field of each subject
    fields{i} = m.node_data{get_field_idx(m, field_name, 'node')}.data;
end
%% Calculate and plot averages
% Calculate
fields = cell2mat(fields);
avg_field = mean(fields, 2);
std_field = std(fields, 0, 2);
% Plot
m.node_data = {}; %cleanup fields
m.node_data{1}.data = avg_field; % add average field
m.node_data{1}.name = [field_name '_avg'];
m.node_data{2}.data = std_field; % add std field
m.node_data{2}.name = [field_name '_std'];

% show surfaces with fields
h5=mesh_show_surfaceycc(m, 'field_idx', [field_name '_avg'],'scaleLimits',[-0.1 0.1]);
title('E_normal_avg Ring')
h6=mesh_show_surfaceycc(m, 'field_idx', [field_name '_std'],'scaleLimits',[-0.1 0.1]);
title('E_normal_std Ring')
savefig(h5,'E_normal_avg Ring')
savefig(h6,'E_normal_std Ring')
save('E_normal_avg_Ring','m')

%%%%%%%%%FILED ENORM
field_name = 'E_norm';
fields = {};
for i = 1:length(subjects)
    sub = subjects{i};
    % load mesh with results transformed to fsaverage space
    m = mesh_load_gmsh4(fullfile(pwd, sub, results_folder, [sub fsavg_msh_name]));
    % Save the field of each subject
    fields{i} = m.node_data{get_field_idx(m, field_name, 'node')}.data;
end
%% Calculate and plot averages
% Calculate
fields = cell2mat(fields);
avg_field = mean(fields, 2);
std_field = std(fields, 0, 2);
% Plot
m.node_data = {}; %cleanup fields
m.node_data{1}.data = avg_field; % add average field
m.node_data{1}.name = [field_name '_avg'];
m.node_data{2}.data = std_field; % add std field
m.node_data{2}.name = [field_name '_std'];

% show surfaces with fields
h7=mesh_show_surfaceycc(m, 'field_idx', [field_name '_avg'],'scaleLimits',[-0.19 0.19]);
title('E_norm_avg Ring')
h8=mesh_show_surfaceycc(m, 'field_idx', [field_name '_std'],'scaleLimits',[-0.1 0.1]);
title('E_norm_std Ring')
savefig(h7,'E_norm_avg Ring')
savefig(h8,'E_norm_std Ring')
save('E_norm_avg_Ring','m')

