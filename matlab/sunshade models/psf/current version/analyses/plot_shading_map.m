%% plot_shading_map.m
%  Reads the irradiance factor NC file and plots a shaded globe of the Earth
%  for a specified timestamp.
%
%  Usage: Run this script directly in MATLAB.
%  The NC file path can be changed in the USER SETTINGS section below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  USER SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nc_file    = fullfile(config_paths().nc_exports, 'v.26.03.19-11.08 irradiance factors.nc');
time_index = 1;   % Which time step to plot (1 = first day)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  READ NC FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lat  = ncread(nc_file, 'lat');          % 192x1, degrees north
lon  = ncread(nc_file, 'lon');          % 288x1, degrees east
df   = ncread(nc_file, 'df');           % 288x192x365 (lon x lat x time)

% Extract the time slice and transpose to lat x lon
df_slice = squeeze(df(:, :, time_index))';   % 192x288 (lat x lon)

% Shift lon from 0-360 to -180-180
lon_shifted = lon;
lon_shifted(lon > 180) = lon(lon > 180) - 360;
[lon_sorted, sort_idx] = sort(lon_shifted);
df_sorted = df_slice(:, sort_idx);

% Read the date string
date_chars = ncread(nc_file, 'date');
date_str   = date_chars(:, time_index)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  BUILD GLOBE GEOMETRY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert lat/lon grid to radians
[lon_grid, lat_grid] = meshgrid(lon_sorted, lat);
lon_rad = deg2rad(lon_grid);
lat_rad = deg2rad(lat_grid);

% Unit sphere XYZ
X = cos(lat_rad) .* cos(lon_rad);
Y = cos(lat_rad) .* sin(lon_rad);
Z = sin(lat_rad);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color', 'k', 'Position', [100 100 900 900]);
clf;
ax = axes('Color', 'k');

%% --- Irradiance factor on the globe surface ---
s = surf(X, Y, Z, df_sorted, 'EdgeColor', 'none');
shading interp;

colormap(ax, interp1([0 1], [0.2 0.4 1; 1 0.5 0], linspace(0,1,256)));
clim([min(df_sorted(:)) max(df_sorted(:))]);
cb = colorbar;
cb.Label.String   = 'Irradiance factor (1 = full sun, 0 = fully shaded)';
cb.Label.FontSize = 11;
cb.Color          = 'w';

hold on;

%% --- Coastline outline ---
coast_url  = 'https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_110m_coastline.geojson';
coast_file = fullfile(tempdir, 'coastline.geojson');
if ~isfile(coast_file)
    websave(coast_file, coast_url);
end
coast_data = jsondecode(fileread(coast_file));
r = 1.001;   % Slightly above the sphere surface so lines are visible
for k = 1:numel(coast_data.features)
    coords    = coast_data.features(k).geometry.coordinates;
    clon_rad  = deg2rad(coords(:,1));
    clat_rad  = deg2rad(coords(:,2));
    cx = r * cos(clat_rad) .* cos(clon_rad);
    cy = r * cos(clat_rad) .* sin(clon_rad);
    cz = r * sin(clat_rad);
    plot3(cx, cy, cz, 'k-', 'LineWidth', 0.6);
end

%% --- Formatting ---
axis equal off;
view(0, 20);   % Viewing angle: [azimuth, elevation] — adjust to taste
lighting phong;
camlight('headlight');

title(sprintf('Planetary Sunshade — Irradiance Factor\nDate: %s   (time index %d)', ...
              date_str, time_index), ...
      'FontSize', 14, 'Color', 'w');

%% --- Summary statistics: area-weighted mean from central column only ---
central_col = ceil(size(df_sorted, 2) / 2);
weights     = cosd(lat)';
weights     = weights / sum(weights);
mean_df     = sum(df_sorted(:, central_col) .* weights);
min_df      = min(df_sorted(:), [], 'omitnan');
frac_shaded = mean(df_sorted(:) < max(df_sorted(:)), 'omitnan') * 100;

annotation('textbox', [0.01 0.01 0.30 0.10], ...
    'String', sprintf('Mean irradiance factor: %.8f\nMin irradiance factor:  %.8f\nFraction of Earth shaded: %.1f%%', ...
                      mean_df, min_df, frac_shaded), ...
    'EdgeColor', 'w', 'BackgroundColor', 'k', 'Color', 'w', ...
    'FontSize', 9, 'FitBoxToText', 'on');

fprintf('\n--- Time index %d  (%s) ---\n', time_index, date_str);
fprintf('Mean irradiance factor : %.8f\n', mean_df);
fprintf('Min  irradiance factor : %.8f\n', min_df);
fprintf('Fraction of globe < 1  : %.1f%%\n', frac_shaded);