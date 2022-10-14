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
erps_age_diff_tonelength_fz = [];
erps_age_diff_tonelength_pz = [];
erps_age_diff_oddball_fz = [];
erps_age_diff_oddball_pz = [];
erps_vnorm_age_diff_tonelength_fz = [];
erps_vnorm_age_diff_tonelength_pz = [];
erps_vnorm_age_diff_oddball_fz = [];
erps_vnorm_age_diff_oddball_pz = [];
erps_maxscal_age_diff_tonelength_fz = [];
erps_maxscal_age_diff_tonelength_pz = [];
erps_maxscal_age_diff_oddball_fz = [];
erps_maxscal_age_diff_oddball_pz = [];
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

    % Difference waves of minmax-scaled ERPs
    erps_maxscal_age_diff_tonelength_fz(a, :) = erp_fz_long / max(abs(erp_fz_long)) - erp_fz_short / max(abs(erp_fz_short));
    erps_maxscal_age_diff_tonelength_pz(a, :) = erp_pz_long / max(abs(erp_pz_long)) - erp_pz_short / max(abs(erp_pz_short));
    erps_maxscal_age_diff_oddball_fz(a, :)    = erp_fz_dev  / max(abs(erp_fz_dev))  - erp_fz_std   / max(abs(erp_fz_std));
    erps_maxscal_age_diff_oddball_pz(a, :)    = erp_pz_dev  / max(abs(erp_pz_dev))  - erp_pz_std   / max(abs(erp_pz_std));

    % Difference waves of vector-normalized ERPs
    erps_vnorm_age_diff_tonelength_fz(a, :) = erp_fz_long / sqrt(sum(erp_fz_long.^2)) - erp_fz_short / sqrt(sum(erp_fz_short.^2));
    erps_vnorm_age_diff_tonelength_pz(a, :) = erp_pz_long / sqrt(sum(erp_pz_long.^2)) - erp_pz_short / sqrt(sum(erp_pz_short.^2));
    erps_vnorm_age_diff_oddball_fz(a, :)    = erp_fz_dev / sqrt(sum(erp_fz_dev.^2)) - erp_fz_std / sqrt(sum(erp_fz_std.^2));
    erps_vnorm_age_diff_oddball_pz(a, :)    = erp_pz_dev / sqrt(sum(erp_pz_dev.^2)) - erp_pz_std / sqrt(sum(erp_pz_std.^2));

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



figure()
subplot(2, 2, 1)
pd = erps_vnorm_age_diff_tonelength_fz;
contourf(erp_times, agebins, pd, 40, 'linecolor','none')
colormap('jet')
set(gca, 'clim', [-max(abs(pd(:))), max(abs(pd(:)))])
colorbar;
xlabel('ms')
ylabel('age')
title(['vnorm - long - short at Fz'])

subplot(2, 2, 2)
pd = erps_vnorm_age_diff_oddball_fz;
contourf(erp_times, agebins, pd, 40, 'linecolor','none')
colormap('jet')
set(gca, 'clim', [-max(abs(pd(:))), max(abs(pd(:)))])
colorbar;
xlabel('ms')
ylabel('age')
title(['vnorm - deviant - standard at Fz'])

subplot(2, 2, 3)
pd = erps_vnorm_age_diff_tonelength_pz;
contourf(erp_times, agebins, pd, 40, 'linecolor','none')
colormap('jet')
set(gca, 'clim', [-max(abs(pd(:))), max(abs(pd(:)))])
colorbar;
xlabel('ms')
ylabel('age')
title(['vnorm - long - short at Pz'])

subplot(2, 2, 4)
pd = erps_vnorm_age_diff_oddball_pz;
contourf(erp_times, agebins, pd, 40, 'linecolor','none')
colormap('jet')
set(gca, 'clim', [-max(abs(pd(:))), max(abs(pd(:)))])
colorbar;
xlabel('ms')
ylabel('age')
title(['vnorm - deviant - standard at Pz'])



figure()
subplot(2, 2, 1)
pd = erps_maxscal_age_diff_tonelength_fz;
contourf(erp_times, agebins, pd, 40, 'linecolor','none')
colormap('jet')
set(gca, 'clim', [-max(abs(pd(:))), max(abs(pd(:)))])
colorbar;
xlabel('ms')
ylabel('age')
title(['maxscal - long - short at Fz'])

subplot(2, 2, 2)
pd = erps_maxscal_age_diff_oddball_fz;
contourf(erp_times, agebins, pd, 40, 'linecolor','none')
colormap('jet')
set(gca, 'clim', [-max(abs(pd(:))), max(abs(pd(:)))])
colorbar;
xlabel('ms')
ylabel('age')
title(['maxscal - deviant - standard at Fz'])

subplot(2, 2, 3)
pd = erps_maxscal_age_diff_tonelength_pz;
contourf(erp_times, agebins, pd, 40, 'linecolor','none')
colormap('jet')
set(gca, 'clim', [-max(abs(pd(:))), max(abs(pd(:)))])
colorbar;
xlabel('ms')
ylabel('age')
title(['maxscal - long - short at Pz'])

subplot(2, 2, 4)
pd = erps_maxscal_age_diff_oddball_pz;
contourf(erp_times, agebins, pd, 40, 'linecolor','none')
colormap('jet')
set(gca, 'clim', [-max(abs(pd(:))), max(abs(pd(:)))])
colorbar;
xlabel('ms')
ylabel('age')
title(['maxscal - deviant - standard at Pz'])

