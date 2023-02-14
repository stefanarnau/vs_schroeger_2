% Residuals...
clear all;

% Path vars
PATH_PREPROCESSED = '/mnt/data_dump/schroeger2/preprocessed/';
PATH_OUT = '/mnt/data_dump/schroeger2/results/';

% Get list of files
fl = dir(fullfile(PATH_PREPROCESSED, '*.mat'));

% Read list of valid subjects
valid_subjects = readmatrix('list_valid_vp.csv');

% Compare eeg to valid-list
eeg_in_valid = [];
eeg_not_in_valid = [];
eeg_all = [];
for f = 1 : numel(fl)
    id = str2num(fl(f).name(4 : 6));
    if ismember(id, valid_subjects)
        eeg_in_valid(end + 1) = id;
    else
        eeg_not_in_valid(end + 1) = id;
    end
    eeg_all(end + 1) = id;
end

% Compare valid-list to eeg data
valid_in_eeg = [];
valid_not_in_eeg = [];
for f = 1 : length(valid_subjects)
    id = valid_subjects(f);
    if ismember(id, eeg_all)
        valid_in_eeg(end + 1) = id;
    else
        valid_not_in_eeg(end + 1) = id;
    end
end

% Save valid in EEG
writematrix(valid_in_eeg, 'valide_in_eeg_data.csv');

% Dummy load for dimensionality
load([PATH_PREPROCESSED, fl(1).name]);

% Get time
erp_times = VP_Header.erp_time;

% Need space!
ages = zeros(length(valid_in_eeg), 1);
ids = [];
erp_data = zeros(length(valid_in_eeg), 4, length(VP_Header.chanlocs), length(VP_Header.erp_time));

% Iterate datasets
idx_nan_data = [];
ids_nan_data = [];
sum_nan_data = [];
counter = 0;
for f = 1 : numel(fl)

    % Talk
    fprintf('\nread dataset %i/%i...\n', f, length(valid_in_eeg));

    % Load data
    load([PATH_PREPROCESSED, fl(f).name]);

    % Get id
    id = str2num(fl(f).name(4 : 6));

    % Check id
    if ~ismember(id, valid_in_eeg)
        continue;
    end

    counter = counter + 1;

    % Get relevant data
    erp_data(counter, :, :, :) = ERP_Dat.data;
    ages(counter, 1) = VP_Header.Age;
    ids(counter, 1) = str2num(VP_Header.ID(4 : 6));

end

% Prune in time
idx_time = erp_times >= -200 & erp_times <= 1000;
erp_times = erp_times(idx_time);
erp_data = erp_data(:, :, :, idx_time);

% Remove mastoid channels
erp_data(:, :, [10, 21], :) = [];
VP_Header.chanlocs([10, 21]) = [];

% New electrode order
new_order =    [1, 29, 30,...          % FP1 FPz Fp2 
                3, 4, 31, 27, 28,...   % F7 F3 Fz F4 F8
                2, 5, 26,...           % FC3 FCz FC4 
                7, 8, 32, 23, 24,...   % T7 C3 Cz C4 T8
                6, 9, 25,...           % CP3 CPz CP4
                11, 12, 13, 19, 20,... % P7 P3 Pz P4 P8
                14, 22, 18,...         % PO3 POz PO4
                15, 16, 17];           % O1 Oz O2

% Adjust new order to removed mastoids
new_order(new_order > 21) = new_order(new_order > 21) - 1;
new_order(new_order > 10) = new_order(new_order > 10) - 1;

% Re-order channels
erp_data = erp_data(:, :, new_order, :);
chanlocs = VP_Header.chanlocs;
for ch = 1 : numel(VP_Header.chanlocs)
    chanlocs(ch) = VP_Header.chanlocs(new_order(ch));
end

% Save erp data [id x condition x channel x time]
save([PATH_OUT, 'erp_data.mat'], 'erp_data', 'ages', 'erp_times');

% Conditions:
% 1: std short
% 2: std long
% 3: dev short
% 4: dev long

% Get min and max age
agelimits = [floor(min(ages)), floor(max(ages))];

% Create agebins
agebins = agelimits(1) : agelimits(end);

% Loop ages in 1 year bins
erps_age_diff_tonelength_fz = [];
erps_age_diff_tonelength_pz = [];
erps_age_diff_oddball_fz = [];
erps_age_diff_oddball_pz = [];
n_bin = [];
for a = 1 : length(agebins)

    % Get age idx
    agesmearing = 1;
    idx_age = ages >= agebins(a) & ages < agebins(a) + agesmearing;

    % Get number of datasets in bin
    n_bin(a) = sum(idx_age);

    % Get age-bin erps
    age_erps = squeeze(mean(erp_data(idx_age, :, :, :), 1));

    % Get age-bin erps at certain electrodes
    age_erps_pz = squeeze(age_erps(:, find(strcmpi({chanlocs.labels}, 'Pz')), :));
    age_erps_fz = squeeze(age_erps(:, find(strcmpi({chanlocs.labels}, 'Fz')), :));

    % Get some ERPs
    erp_fz_long  = mean(age_erps_fz([2, 4], :), 1);
    erp_fz_short = mean(age_erps_fz([1, 3], :), 1);
    erp_pz_long  = mean(age_erps_pz([2, 4], :), 1);
    erp_pz_short = mean(age_erps_pz([1, 3], :), 1);
    erp_fz_std   = mean(age_erps_fz([1, 2], :), 1);
    erp_fz_dev   = mean(age_erps_fz([3, 4], :), 1);
    erp_pz_std   = mean(age_erps_pz([1, 2], :), 1);
    erp_pz_dev   = mean(age_erps_pz([3, 4], :), 1);

    % Difference waves
    erps_age_diff_tonelength_fz(a, :) = erp_fz_long - erp_fz_short;
    erps_age_diff_tonelength_pz(a, :) = erp_pz_long - erp_pz_short;
    erps_age_diff_oddball_fz(a, :) = erp_fz_dev - erp_fz_std;
    erps_age_diff_oddball_pz(a, :) = erp_pz_dev - erp_pz_std;

end

figure()
subplot(2, 2, 1)
pd = erps_age_diff_tonelength_fz;
contourf(erp_times, agebins, pd, 40, 'linecolor','none')
colormap('jet')
set(gca, 'clim', [-max(abs(pd(:))), max(abs(pd(:)))])
colorbar;
xlabel('ms')
ylabel('age')
title(['long - short at Fz'])

subplot(2, 2, 2)
pd = erps_age_diff_oddball_fz;
contourf(erp_times, agebins, pd, 40, 'linecolor','none')
colormap('jet')
set(gca, 'clim', [-max(abs(pd(:))), max(abs(pd(:)))])
colorbar;
xlabel('ms')
ylabel('age')
title(['deviant - standard at Fz'])

subplot(2, 2, 3)
pd = erps_age_diff_tonelength_pz;
contourf(erp_times, agebins, pd, 40, 'linecolor','none')
colormap('jet')
set(gca, 'clim', [-max(abs(pd(:))), max(abs(pd(:)))])
colorbar;
xlabel('ms')
ylabel('age')
title(['long - short at Pz'])

subplot(2, 2, 4)
pd = erps_age_diff_oddball_pz;
contourf(erp_times, agebins, pd, 40, 'linecolor','none')
colormap('jet')
set(gca, 'clim', [-max(abs(pd(:))), max(abs(pd(:)))])
colorbar;
xlabel('ms')
ylabel('age')
title(['deviant - standard at Pz'])

% Init ft
PATH_FIELDTRIP = '/home/plkn/fieldtrip-master/';
addpath(PATH_FIELDTRIP);
ft_defaults;

% Read example header for channel labels
hdr = ft_read_header('/mnt/data_dump/schroeger2/DKA001_Schroeger_ICA.set')

% Init EEGLab
PATH_EEGLAB = '/home/plkn/eeglab2022.1/';
addpath(PATH_EEGLAB);
eeglab;

% Get difference waves (dims: subject x condition x channel x time)
erp_std = squeeze(mean(erp_data(:, [1, 2], :, :), 2));
erp_dev = squeeze(mean(erp_data(:, [3, 4], :, :), 2));

% Difference dims: subject x channel x time
erp_diff = erp_dev - erp_std;

% Bild elec struct
elec = struct();
for ch = 1 : length(chanlocs)
    elec.label{ch} = chanlocs(ch).labels;
    elec.elecpos(ch, :) = [chanlocs(ch).X, chanlocs(ch).Y, chanlocs(ch).Z];
    elec.chanpos(ch, :) = [chanlocs(ch).X, chanlocs(ch).Y, chanlocs(ch).Z];
end

% Build struct for ft
for s = 1 : length(ages)
    d = [];
    d.dimord = 'chan_time';
    d.label = elec.label;
    d.time = erp_times;
    d.avg = squeeze(erp_diff(s, :, :));
    D{s} = d;
end

cfg=[];
cfg.keepindividual = 'yes';
GA = ft_timelockgrandaverage(cfg, D{1, :});

% Prepare layout
cfg      = [];
cfg.elec = elec;
cfg.rotate = 90;
layout = ft_prepare_layout(cfg);

% Define neighbours
cfg                 = [];
cfg.layout          = layout;
cfg.feedback        = 'no';
cfg.method          = 'triangulation'; 
cfg.neighbours      = ft_prepare_neighbours(cfg, GA);
neighbours = cfg.neighbours;

% The design
design = ages;

% Set config
cfg = [];

% One sided test!
cfg.tail             = 0; 

% Made changes to ft_statfun_correlationT.m. Line 78 design->design(cfg.ivar,:)
cfg.statistic        = 'ft_statfun_correlationT'; 

cfg.alpha            = 0.025;
cfg.neighbours       = neighbours;
cfg.minnbchan        = 2;
cfg.method           = 'montecarlo';
cfg.correctm         = 'cluster';
cfg.clustertail      = 0;
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.numrandomization = 1000;
cfg.computecritval   = 'yes'; 
cfg.ivar             = 1;
cfg.design           = design;

% The test
[stat] = ft_timelockstatistics(cfg, GA);  

% Indices of the significant cluster
idx = stat.posclusterslabelmat == 1;

% Identify significant channels and time points
chans_sig = find(sum(idx, 2));
times_sig = find(sum(idx, 1));

aa = bb;

% Plot a topo for significat time window
markercolor = 'k';
cmap = 'jet';
clim = [-0.5, 0.5];
figure('Visible', 'off'); clf;
pd = mean(stat.rho(:, times_sig), 2);
topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {chans_sig, 'p', markercolor, 14, 1});
colormap(cmap);
caxis(clim);
saveas(gcf, [PATH_OUT 'rho_timewin_topo.png']);

% Save correlation matrix and cluster overlay
dlmwrite([PATH, 'rho_chan_time.csv'], stat.rho);
dlmwrite([PATH, 'cluster_contour_chan_time.csv'], idx);

% Save Line plot averaged across significant channels + axis
dlmwrite([PATH, 'rho_time.csv'], mean(stat.rho(chans_sig, :), 1));
dlmwrite([PATH, 'timevec.csv'], [1 : 1000]);