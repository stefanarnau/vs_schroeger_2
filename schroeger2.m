% Residuals...
clear all;

% Path vars
PATH_PREPROCESSED = '/mnt/data_dump/schroeger2/preprocessed/';

% Get list of files
fl = dir(fullfile(PATH_PREPROCESSED, '*.mat'));

% Dummy load for dimensionality
load([PATH_PREPROCESSED, fl(1).name]);

% Get time
erp_times = VP_Header.erp_time;

% Get chanlocs
chanlocs = VP_Header.chanlocs;

% Need space!
ages = zeros(numel(fl), 1);
erp_data = zeros(numel(fl), 4, length(VP_Header.chanlocs), length(VP_Header.erp_time));

% Iterate datasets
idx_nan_data = [];
sum_nan_data = [];
for f = 1 : numel(fl)

    % Talk
    fprintf('\nread dataset %i/%i...\n', f, numel(fl));

    % Load data
    load([PATH_PREPROCESSED, fl(f).name]);

    % Get relevant data
    erp_data(f, :, :, :) = ERP_Dat.data;
    ages(f, 1) = VP_Header.Age;

    % Check for missing data
    if sum(isnan(ERP_Dat.data), [1, 2, 3])
        idx_nan_data(end + 1) = f;
        sum_nan_data(end + 1) = sum(isnan(ERP_Dat.data), [1, 2, 3]);
    end

end

% Delete missing data datasets
ages(idx_nan_data) = [];
erp_data(idx_nan_data, :, :, :) = [];

% Prune in time
idx_time = erp_times >= -200 & erp_times <= 1000;
erp_times = erp_times(idx_time);
erp_data = erp_data(:, :, :, idx_time);

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
erps_age_diff_tonelength = [];
erps_age_diff_oddball = [];
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
    age_erps_pz = squeeze(age_erps(:, 13, :));
    age_erps_fz = squeeze(age_erps(:, 31, :));

    % Tone length difference wave
    erps_age_diff_tonelength(a, :) = mean(age_erps_pz([2, 4], :), 1) - mean(age_erps_pz([1, 3], :), 1);

    % Oddball difference wave
    erps_age_diff_oddball(a, :) = mean(age_erps_fz([3, 4], :), 1) - mean(age_erps_fz([1, 2], :), 1);

end

figure()
subplot(1, 2, 1)
pd = erps_age_diff_tonelength;
contourf(erp_times, agebins, pd, 40, 'linecolor','none')
colormap('jet')
set(gca, 'clim', [-max(abs(pd(:))), max(abs(pd(:)))])
colorbar;
xlabel('ms')
ylabel('age')
title(['long - short at Pz'])

subplot(1, 2, 2)
pd = erps_age_diff_oddball;
contourf(erp_times, agebins, pd, 40, 'linecolor','none')
colormap('jet')
set(gca, 'clim', [-max(abs(pd(:))), max(abs(pd(:)))])
colorbar;
xlabel('ms')
ylabel('age')
title(['deviant - standard at Fz'])











