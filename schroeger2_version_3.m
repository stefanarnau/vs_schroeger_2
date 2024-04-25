% Residuals...
clear all;

% Path vars
PATH_PREPROCESSED = '/mnt/data_dump/schroeger2/preprocessed/';
PATH_OUT = '/mnt/data_dump/schroeger2/results/';

% Get list of files
fl = dir(fullfile(PATH_PREPROCESSED, '*.mat'));

% Read list of valid subjects
valid_subjects = readmatrix('Schröger-VP_24-04-09.csv');

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
erps_age_diff_frontal_oddball_short = [];
erps_age_diff_frontal_oddball_long  = [];
erps_age_diff_parietal_oddball_short = [];
erps_age_diff_parietal_oddball_long  = [];

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
    channel_idx_frontal  = [1, 2, 3, 5, 6, 7, 9, 10, 11]; 
    channel_idx_parietal = [21, 22, 23, 25, 26, 27];
    age_erps_frontal  = squeeze(mean(age_erps(:, channel_idx_frontal,  :), 2));
    age_erps_parietal = squeeze(mean(age_erps(:, channel_idx_parietal, :), 2));

    % Difference waves
    erps_age_diff_frontal_oddball_short(a, :) = age_erps_frontal(3, :) - age_erps_frontal(1, :);
    erps_age_diff_frontal_oddball_long(a, :)  = age_erps_frontal(4, :) - age_erps_frontal(2, :);

    erps_age_diff_parietal_oddball_short(a, :) = age_erps_parietal(3, :) - age_erps_parietal(1, :);
    erps_age_diff_parietal_oddball_long(a, :)  = age_erps_parietal(4, :) - age_erps_parietal(2, :);
    
end

% Save diff-erps for plotting
dlmwrite([PATH_OUT, 'erps_age_diff_frontal_oddball_short.csv'], erps_age_diff_frontal_oddball_short);
dlmwrite([PATH_OUT, 'erps_age_diff_frontal_oddball_long.csv'], erps_age_diff_frontal_oddball_long);
dlmwrite([PATH_OUT, 'erps_age_diff_parietal_oddball_short.csv'], erps_age_diff_parietal_oddball_short);
dlmwrite([PATH_OUT, 'erps_age_diff_parietal_oddball_long.csv'], erps_age_diff_parietal_oddball_long);

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

% Get difference waves deviant - standard for short and for long (dims: subject x channel x time)
erp_diff_srt = squeeze(erp_data(:, 3, :, :)) - squeeze(erp_data(:, 1, :, :));
erp_diff_lng = squeeze(erp_data(:, 4, :, :)) - squeeze(erp_data(:, 2, :, :));

% Build elec struct
elec = struct();
for ch = 1 : length(chanlocs)
    elec.label{ch} = chanlocs(ch).labels;
    elec.elecpos(ch, :) = [chanlocs(ch).X, chanlocs(ch).Y, chanlocs(ch).Z];
    elec.chanpos(ch, :) = [chanlocs(ch).X, chanlocs(ch).Y, chanlocs(ch).Z];
end

% Build GA dev minus standard for short and for long
for s = 1 : length(ages)
    d = [];
    d.dimord = 'chan_time';
    d.label = elec.label;
    d.time = erp_times;
    d.avg = squeeze(erp_diff_srt(s, :, :));
    D{s} = d;
end
cfg=[];
cfg.keepindividual = 'yes';
GA_diff_srt = ft_timelockgrandaverage(cfg, D{1, :});
for s = 1 : length(ages)
    d = [];
    d.dimord = 'chan_time';
    d.label = elec.label;
    d.time = erp_times;
    d.avg = squeeze(erp_diff_lng(s, :, :));
    D{s} = d;
end
cfg=[];
cfg.keepindividual = 'yes';
GA_diff_lng = ft_timelockgrandaverage(cfg, D{1, :});

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
cfg.neighbours      = ft_prepare_neighbours(cfg, GA_diff_srt);
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
cfg.clusteralpha     = 0.05; % 0.05 for main effects good!
cfg.clusterstatistic = 'maxsum';
cfg.numrandomization = 1000;
cfg.computecritval   = 'yes'; 
cfg.ivar             = 1;
cfg.design           = design;

% The tests
[stat_srt] = ft_timelockstatistics(cfg, GA_diff_srt);  
[stat_lng] = ft_timelockstatistics(cfg, GA_diff_lng);  

% Save correlation matrix
dlmwrite([PATH_OUT, 'rho_chan_time_srt.csv'], stat_srt.rho);
dlmwrite([PATH_OUT, 'rho_chan_time_lng.csv'], stat_lng.rho);

% Short: Cluster pos 1
idx = stat_srt.posclusterslabelmat == 1;
chans_sig = find(sum(idx, 2));
times_sig = find(sum(idx, 1));
markercolor = 'k';
cmap = 'jet';
clim = [-0.5, 0.5];
figure('Visible', 'off'); clf;
pd = mean(stat_srt.rho(:, times_sig), 2);
topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {chans_sig, '.', markercolor, 14, 1});
colormap(cmap);
caxis(clim);
saveas(gcf, [PATH_OUT 'topo_cluster_short_pos_1.png']);
dlmwrite([PATH_OUT, 'contour_cluster_short_pos_1.csv'], idx);

% Short: Cluster pos 2
idx = stat_srt.posclusterslabelmat == 2;
chans_sig = find(sum(idx, 2));
times_sig = find(sum(idx, 1));
markercolor = 'k';
cmap = 'jet';
clim = [-0.5, 0.5];
figure('Visible', 'off'); clf;
pd = mean(stat_srt.rho(:, times_sig), 2);
topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {chans_sig, '.', markercolor, 14, 1});
colormap(cmap);
caxis(clim);
saveas(gcf, [PATH_OUT 'topo_cluster_short_pos_2.png']);
dlmwrite([PATH_OUT, 'contour_cluster_short_pos_2.csv'], idx);

% Short: Cluster neg 1
idx = stat_srt.negclusterslabelmat == 1;
chans_sig = find(sum(idx, 2));
times_sig = find(sum(idx, 1));
markercolor = 'k';
cmap = 'jet';
clim = [-0.5, 0.5];
figure('Visible', 'off'); clf;
pd = mean(stat_srt.rho(:, times_sig), 2);
topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {chans_sig, '.', markercolor, 14, 1});
colormap(cmap);
caxis(clim);
saveas(gcf, [PATH_OUT 'topo_cluster_short_neg_1.png']);
dlmwrite([PATH_OUT, 'contour_cluster_short_neg_1.csv'], idx);

% Short: Cluster neg 2
idx = stat_srt.negclusterslabelmat == 2;
chans_sig = find(sum(idx, 2));
times_sig = find(sum(idx, 1));
markercolor = 'k';
cmap = 'jet';
clim = [-0.5, 0.5];
figure('Visible', 'off'); clf;
pd = mean(stat_srt.rho(:, times_sig), 2);
topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {chans_sig, '.', markercolor, 14, 1});
colormap(cmap);
caxis(clim);
saveas(gcf, [PATH_OUT 'topo_cluster_short_neg_2.png']);
dlmwrite([PATH_OUT, 'contour_cluster_short_neg_2.csv'], idx);

% long: Cluster pos 1
idx = stat_lng.posclusterslabelmat == 1;
chans_sig = find(sum(idx, 2));
times_sig = find(sum(idx, 1));
markercolor = 'k';
cmap = 'jet';
clim = [-0.5, 0.5];
figure('Visible', 'off'); clf;
pd = mean(stat_lng.rho(:, times_sig), 2);
topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {chans_sig, '.', markercolor, 14, 1});
colormap(cmap);
caxis(clim);
saveas(gcf, [PATH_OUT 'topo_cluster_long_pos_1.png']);
dlmwrite([PATH_OUT, 'contour_cluster_long_pos_1.csv'], idx);

% Long: Cluster neg 1
idx = stat_lng.negclusterslabelmat == 1;
chans_sig = find(sum(idx, 2));
times_sig = find(sum(idx, 1));
markercolor = 'k';
cmap = 'jet';
clim = [-0.5, 0.5];
figure('Visible', 'off'); clf;
pd = mean(stat_lng.rho(:, times_sig), 2);
topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {chans_sig, '.', markercolor, 14, 1});
colormap(cmap);
caxis(clim);
saveas(gcf, [PATH_OUT 'topo_cluster_long_neg_1.png']);
dlmwrite([PATH_OUT, 'contour_cluster_long_neg_1.csv'], idx);

% Long: Cluster neg 2
idx = stat_lng.negclusterslabelmat == 2;
chans_sig = find(sum(idx, 2));
times_sig = find(sum(idx, 1));
markercolor = 'k';
cmap = 'jet';
clim = [-0.5, 0.5];
figure('Visible', 'off'); clf;
pd = mean(stat_lng.rho(:, times_sig), 2);
topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {chans_sig, '.', markercolor, 14, 1});
colormap(cmap);
caxis(clim);
saveas(gcf, [PATH_OUT 'topo_cluster_long_neg_2.png']);
dlmwrite([PATH_OUT, 'contour_cluster_long_neg_2.csv'], idx);




















% pos cluster 2
idx = stat_dev_vs_std.posclusterslabelmat == 2;

% Identify significant channels and time points
chans_sig = find(sum(idx, 2));
times_sig = find(sum(idx, 1));

% Plot a topo for significat time window
markercolor = 'k';
cmap = 'jet';
clim = [-0.5, 0.5];
figure('Visible', 'off'); clf;
pd = mean(stat_dev_vs_std.rho(:, times_sig), 2);
topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {chans_sig, '.', markercolor, 14, 1});
colormap(cmap);
caxis(clim);
saveas(gcf, [PATH_OUT 'rho_timewin_topo_dev_vs_std_pos_2.png']);

% Save contour
dlmwrite([PATH_OUT, 'cluster_contour_dev_vs_std_pos_2.csv'], idx);

% neg cluster 1
idx = stat_dev_vs_std.negclusterslabelmat == 1;

% Identify significant channels and time points
chans_sig = find(sum(idx, 2));
times_sig = find(sum(idx, 1));

% Plot a topo for significat time window
markercolor = 'k';
cmap = 'jet';
clim = [-0.5, 0.5];
figure('Visible', 'off'); clf;
pd = mean(stat_dev_vs_std.rho(:, times_sig), 2);
topoplot(pd, chanlocs, 'plotrad', 0.7, 'intrad', 0.7, 'intsquare', 'on', 'conv', 'off', 'electrodes', 'off', 'emarker2', {chans_sig, '.', markercolor, 14, 1});
colormap(cmap);
caxis(clim);
saveas(gcf, [PATH_OUT 'rho_timewin_topo_dev_vs_std_neg_1.png']);

% Save contour
dlmwrite([PATH_OUT, 'cluster_contour_dev_vs_std_neg_1.csv'], idx);

% Correlation lineplots =================================================================================

% Init matrices
erp_Fz_std_short = zeros(size(erp_data, 1), size(erp_data, 4));
erp_Fz_std_long  = zeros(size(erp_data, 1), size(erp_data, 4));
erp_Fz_dev_short = zeros(size(erp_data, 1), size(erp_data, 4));
erp_Fz_dev_long  = zeros(size(erp_data, 1), size(erp_data, 4));
erp_Pz_std_short = zeros(size(erp_data, 1), size(erp_data, 4));
erp_Pz_std_long  = zeros(size(erp_data, 1), size(erp_data, 4));
erp_Pz_dev_short = zeros(size(erp_data, 1), size(erp_data, 4));
erp_Pz_dev_long  = zeros(size(erp_data, 1), size(erp_data, 4));

erp_Fz_diff_short = zeros(size(erp_data, 1), size(erp_data, 4));
erp_Fz_diff_long  = zeros(size(erp_data, 1), size(erp_data, 4));

% erp_data is: subject x condition x channel x time
for s = 1 : size(erp_data, 1)
    
    % Fill
    erp_Fz_std_short(s, :) = squeeze(erp_data(s, 1, find(strcmpi({chanlocs.labels}, 'FCz')), :));
    erp_Fz_std_long(s, :)  = squeeze(erp_data(s, 2, find(strcmpi({chanlocs.labels}, 'FCz')), :));
    erp_Fz_dev_short(s, :) = squeeze(erp_data(s, 3, find(strcmpi({chanlocs.labels}, 'FCz')), :));
    erp_Fz_dev_long(s, :)  = squeeze(erp_data(s, 4, find(strcmpi({chanlocs.labels}, 'FCz')), :));
    erp_Pz_std_short(s, :) = squeeze(erp_data(s, 1, find(strcmpi({chanlocs.labels}, 'Pz')), :));
    erp_Pz_std_long(s, :)  = squeeze(erp_data(s, 2, find(strcmpi({chanlocs.labels}, 'Pz')), :));
    erp_Pz_dev_short(s, :) = squeeze(erp_data(s, 3, find(strcmpi({chanlocs.labels}, 'Pz')), :));
    erp_Pz_dev_long(s, :)  = squeeze(erp_data(s, 4, find(strcmpi({chanlocs.labels}, 'Pz')), :));

    erp_Fz_diff_short(s, :) = squeeze(erp_data(s, 3, find(strcmpi({chanlocs.labels}, 'Fz')), :)) - squeeze(erp_data(s, 1, find(strcmpi({chanlocs.labels}, 'Fz')), :));
    erp_Fz_diff_long(s, :)  = squeeze(erp_data(s, 4, find(strcmpi({chanlocs.labels}, 'Fz')), :)) - squeeze(erp_data(s, 2, find(strcmpi({chanlocs.labels}, 'Fz')), :));

end

% Correlate
erp_Fz_std_short = corr(erp_Fz_std_short, ages);
erp_Fz_std_long  = corr(erp_Fz_std_long,  ages);
erp_Fz_dev_short = corr(erp_Fz_dev_short, ages);
erp_Fz_dev_long  = corr(erp_Fz_dev_long,  ages);
erp_Pz_std_short = corr(erp_Pz_std_short, ages);
erp_Pz_std_long  = corr(erp_Pz_std_long,  ages);
erp_Pz_dev_short = corr(erp_Pz_dev_short, ages);
erp_Pz_dev_long  = corr(erp_Pz_dev_long,  ages);

erp_Fz_diff_short = corr(erp_Fz_diff_short, ages);
erp_Fz_diff_long  = corr(erp_Fz_diff_long,  ages);

figure()
subplot(1, 2, 1)
plot(erp_times, erp_Fz_std_short)
hold on;
plot(erp_times, erp_Fz_std_long)
plot(erp_times, erp_Fz_dev_short)
plot(erp_times, erp_Fz_dev_long)
legend({'std sh', 'std ln', 'dev sh' 'dev  ln'})
title('FCz')
subplot(1, 2, 2)
plot(erp_times, erp_Pz_std_short)
hold on;
plot(erp_times, erp_Pz_std_long)
plot(erp_times, erp_Pz_dev_short)
plot(erp_times, erp_Pz_dev_long)
legend({'std sh', 'std ln', 'dev sh' 'dev  ln'})
title('Pz')


figure()
plot(erp_times, erp_Fz_diff_short)
hold on;
plot(erp_times, erp_Fz_diff_long)
legend({'sh', 'ln'})
title('FCz')
































































