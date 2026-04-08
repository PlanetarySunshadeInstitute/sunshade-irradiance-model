%% plot_shading_map.m
%  Reads the irradiance factor NC file and plots a shaded globe of the Earth
%  for a specified timestamp.
%
%  The NC file is produced by analysis_General___0___X.m (preconfigured mode).
%  As of the declination-shift fix (April 2026), the 'df' variable encodes
%  geographic latitude correctly: the shadow footprint tracks the solar
%  sub-stellar point across the year rather than being pinned to the equator.
%
%  Verification aid: the title and terminal output include the solar
%  declination delta(day) computed from Spencer (1971).  The sub-stellar
%  latitude is also drawn as a dashed line on the globe.  The shadow trough
%  in the IF map should align with this line.
%
%  Usage: Run this script directly in MATLAB.
%  The NC file path can be changed in the USER SETTINGS section below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  USER SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nc_file    = fullfile(config_paths().nc_exports, 'L1-RX-2026-04-07-001.nc');
time_index = 180;   % Which time step to plot (1 = first day of the year)

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

% Solar declination for this day (Spencer 1971).
% time_index corresponds to day-of-year for a standard annual NC file.
delta_deg  = spencer_declination(time_index, 365);

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

%% --- Sub-stellar latitude ring ---
%  The shadow trough should align with this line after the declination fix.
%  Drawn as a full latitude circle at delta_deg on the sphere surface.
r_ring   = 1.002;
az_ring  = linspace(0, 2*pi, 361);
lat_ring = deg2rad(delta_deg);
cx_ring  = r_ring * cos(lat_ring) .* cos(az_ring);
cy_ring  = r_ring * cos(lat_ring) .* sin(az_ring);
cz_ring  = r_ring * sin(lat_ring) * ones(size(az_ring));
plot3(cx_ring, cy_ring, cz_ring, '--', 'Color', [1 0.85 0], 'LineWidth', 1.2);

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

title(sprintf('Planetary Sunshade — Irradiance Factor\nDate: %s   (day %d)     \delta = %+.2f°  (dashed ring = sub-stellar latitude)', ...
              date_str, time_index, delta_deg), ...
      'FontSize', 13, 'Color', 'w');

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
fprintf('Solar declination      : %+.3f deg  (expected IF trough at lat %+.1f deg)\n', delta_deg, delta_deg);
fprintf('Mean irradiance factor : %.8f\n', mean_df);
fprintf('Min  irradiance factor : %.8f\n', min_df);
fprintf('Fraction of globe < 1  : %.1f%%\n', frac_shaded);