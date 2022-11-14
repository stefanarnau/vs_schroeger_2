% Residuals...
clear all;

% Path vars
PATH_PREPROCESSED = '/mnt/data_dump/schroeger2/preprocessed/';
PATH_OUT = '/mnt/data_dump/schroeger2/results/';

% Get list of files
fl = dir(fullfile(PATH_PREPROCESSED, '*.mat'));

% Dummy load for dimensionality
load([PATH_PREPROCESSED, fl(1).name]);

% Get time
erp_times = VP_Header.erp_time;

% Need space!
ages = zeros(numel(fl), 1);
ids = [];
erp_data = zeros(numel(fl), 4, length(VP_Header.chanlocs), length(VP_Header.erp_time));

% Iterate datasets
idx_nan_data = [];
ids_nan_data = [];
sum_nan_data = [];
for f = 1 : numel(fl)

    % Talk
    fprintf('\nread dataset %i/%i...\n', f, numel(fl));

    % Load data
    load([PATH_PREPROCESSED, fl(f).name]);

    % Get relevant data
    erp_data(f, :, :, :) = ERP_Dat.data;
    ages(f, 1) = VP_Header.Age;
    ids(f, 1) = str2num(VP_Header.ID(4 : 6));

    % Check for missing data
    if sum(isnan(ERP_Dat.data), [1, 2, 3])
        idx_nan_data(end + 1) = f;
        ids_nan_data(end + 1) = str2num(VP_Header.ID(4 : 6));
        sum_nan_data(end + 1) = sum(isnan(ERP_Dat.data), [1, 2, 3]);
    end

end

% Delete missing data datasets
ages(idx_nan_data) = [];
ids(idx_nan_data) = [];
erp_data(idx_nan_data, :, :, :) = [];

% Save ids of in-sample datsets
writematrix(ids, [PATH_OUT, 'ids_of_erp_sample.csv']);

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




