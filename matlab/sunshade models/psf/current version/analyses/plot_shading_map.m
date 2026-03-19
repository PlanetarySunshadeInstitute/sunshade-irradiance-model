%% plot_shading_map.m
%  Reads the irradiance factor NC file and plots a shaded map of the Earth
%  for the first timestamp (time index 1).
%
%  Usage: Run this script directly in MATLAB.
%  The NC file path can be changed in the USER SETTINGS section below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  USER SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nc_file    = '/Users/morgangoodwin/Desktop/PSF/MatLab - Planetary Sunshade Foundation/numerical control/exports/v.26.03.19-06.49 irradiance factors.nc';
time_index = 1;   % Which time step to plot (1 = first day)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  READ NC FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lat  = ncread(nc_file, 'lat');          % 192x1, degrees north
lon  = ncread(nc_file, 'lon');          % 288x1, degrees east
df   = ncread(nc_file, 'df');           % 288x192x365 (lon x lat x time)

% Extract the first time slice and transpose to lat x lon for plotting
df_slice = squeeze(df(:, :, time_index))';   % now 192x288 (lat x lon)

% Read the date string for this time step
date_chars = ncread(nc_file, 'date');        % 5x365 char array
date_str   = date_chars(:, time_index)';     % e.g. '07-04'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color', 'w', 'Position', [100 100 1200 600]);

%% --- Main shading map ---
ax = axes;

% Plot the dimming factor as a filled colour map
% Shift lon from 0-360 to -180-180
lon_shifted = lon;
lon_shifted(lon > 180) = lon(lon > 180) - 360;
[lon_sorted, sort_idx] = sort(lon_shifted);
df_sorted = df_slice(:, sort_idx);
pcolor(lon_sorted, lat, df_sorted);
shading flat;

% Colour scale: 1 = no dimming (full sunlight), 0 = fully shaded
colormap(ax, interp1([0 1], [0.2 0.4 1; 1 0.5 0], linspace(0,1,256)));
clim([min(df_sorted(:)) max(df_sorted(:))]);
cb = colorbar;
cb.Label.String = 'Irradiance factor (1 = full sun, 0 = fully shaded)';
cb.Label.FontSize = 11;

hold on;

%% --- Overlay coastlines ---
% Download coastline data (only needed once)
coast_url  = 'https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_110m_coastline.geojson';
coast_file = fullfile(tempdir, 'coastline.geojson');
if ~isfile(coast_file)
    websave(coast_file, coast_url);
end
coast_data = jsondecode(fileread(coast_file));
for k = 1:numel(coast_data.features)
    coords = coast_data.features(k).geometry.coordinates;
    plot(coords(:,1), coords(:,2), 'k-', 'LineWidth', 0.6);
end

%% --- Formatting ---
axis([-180 180 -90 90]);
xticks(-180:30:180);
yticks(-90:30:90);
xticklabels({'180°W','150°W','120°W','90°W','60°W','30°W','0°', ...
             '30°E','60°E','90°E','120°E','150°E','180°E'});
yticklabels({'90°S','60°S','30°S','0°','30°N','60°N','90°N'});
xlabel('Longitude', 'FontSize', 11);
ylabel('Latitude',  'FontSize', 11);

title(sprintf('Planetary Sunshade — Irradiance Factor\nDate: %s   (time index %d)', ...
              date_str, time_index), 'FontSize', 14);

grid on;
ax.GridColor      = [0.5 0.5 0.5];
ax.GridAlpha      = 0.3;
ax.Layer          = 'top';
ax.FontSize       = 10;
ax.DataAspectRatio = [1 1 1];

%% --- Summary statistics in a text box ---
mean_df = mean(df_sorted(:), 'omitnan');
min_df  = min(df_sorted(:),  [], 'omitnan');
frac_shaded = mean(df_sorted(:) < max(df_sorted(:)), 'omitnan') * 100;

annotation('textbox', [0.01 0.01 0.25 0.10], ...
    'String', sprintf('Mean irradiance factor: %.8f\nMin irradiance factor:  %.8f\nFraction of Earth shaded: %.1f%%', ...
                      mean_df, min_df, frac_shaded), ...
    'EdgeColor', 'k', 'BackgroundColor', 'w', ...
    'FontSize', 9, 'FitBoxToText', 'on');

fprintf('\n--- Time index %d  (%s) ---\n', time_index, date_str);
fprintf('Mean irradiance factor : %.8f\n', mean_df);
fprintf('Min  irradiance factor : %.8f\n', min_df);
fprintf('Fraction of globe < 1  : %.1f%%\n', frac_shaded);