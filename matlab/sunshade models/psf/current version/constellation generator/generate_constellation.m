% generate_constellation.m
%
% Master script for the PSF constellation generator.
% Edit the USER SETTINGS section, then run.
%
% Sequence:
%   1. Load pre-computed Glasgow envelope
%   2. Define and validate constellation parameters
%   3. Place spacecraft
%   4. Visualize — user confirms before export
%   5. Export to excel/psf model/
%
% Authors: Planetary Sunshade Foundation

clear; clc;

% =========================================================================
% USER SETTINGS
% =========================================================================

% --- Paths ---------------------------------------------------------------
envelope_file = fullfile('/Users/morgangoodwin/Desktop/PSF/MatLab/matlab/sunshade models/psf/current version/constellation generator', ...
                          'equilibrium_envelope.mat');
excel_folder  = '/Users/morgangoodwin/Desktop/PSF/MatLab/excel/psf model';

% --- Constellation parameters --------------------------------------------
user_params.pattern                 = 'polar';   % 'uniform' or 'polar'
user_params.N                       = 11003;   % used for uniform mode
user_params.n_planes                = 5;
user_params.plane_spacing_km        = 10000;

user_params.constellation_radius_km = 15000;    % uniform mode only
user_params.footprint_profile       = 'gaussian'; % 'uniform' or 'gaussian'
user_params.footprint_sigma_fraction= 0.1;

user_params.sail_radius_km          = 20;
user_params.min_buffer_km           = 45;

% Polar only (ignored for uniform):
% Specify the number of shades in each polar ellipse explicitly.
user_params.polar_N_north         = 5502;
user_params.polar_N_south         = 5502;
user_params.polar_center_z_km      = 7200;
user_params.polar_radius_y_km      = 3000;    % semi-axis in Y for each pole ellipse
user_params.polar_radius_z_km      = 2000;    % semi-axis in Z for each pole ellipse

% For polar mode, total N is derived from polar_N_north + polar_N_south.
if strcmp(user_params.pattern, 'polar')
    user_params.N = user_params.polar_N_north + user_params.polar_N_south;
end

% =========================================================================
% RUN — no edits needed below this line
% =========================================================================

envelope  = load_envelope(envelope_file);
params    = define_constellation_params(user_params);
[positions, n_violations] = place_spacecraft(params, envelope);

% Pre-compute expected output filename so the stats panel can display it.
params.output_filename = sprintf('constellation_%s_%dcrafts_%s.xlsx', ...
                                  params.pattern, params.N, ...
                                  datestr(now, 'yyyy-mm-dd'));

visualize_constellation(positions, params, envelope, n_violations);

% Pause for user confirmation before writing the Excel file
input_val = input('Constellation looks correct? Export to Excel? (y/n): ', 's');
if ~strcmpi(input_val, 'y')
    fprintf('Export cancelled. Adjust parameters and re-run.\n');
    return;
end

resolved_name = resolve_unique_excel_filename(excel_folder, params.output_filename);
params.output_filename = resolved_name;
filepath = export_kinematics_xlsx(positions, params, excel_folder);

% Update the irradiance model location file to point at the new export.
location_file = fullfile('/Users/morgangoodwin/Desktop/PSF/MatLab/matlab/sunshade models/psf/current version/file locations', ...
                         'location_Heliogyro_Kinematics_Data___E___Sr.m');
update_heliogyro_location_file(location_file, resolved_name);

% Auto-run irradiance analysis and capture its 11x11 grid.
results_grid = analysis_General___0___X;
results_grid = analysis_General___0___X;
% DIAGNOSTIC
disp('Raw grid (row 1 = first model point, row 11 = last):');
disp(results_grid);
if isempty(results_grid) || ~ismatrix(results_grid)
    error('analysis_General___0___X did not return a 2-D irradiance matrix.');
end
fprintf('\n--- Irradiance factor grid (11x11) ---\n');
disp(round(results_grid, 4));

% Replace panel 1 in the existing diagnostic figure with the irradiance grid.
fig = findobj('Type', 'figure', 'Name', 'Constellation Diagnostic');
if isempty(fig)
    fig = gcf;
else
    fig = fig(1);
end
render_irradiance_grid_panel(fig, results_grid);

% Save updated diagnostic figure next to the .xlsx with "_diagnostic" suffix.
[excel_dir, excel_base, ~] = fileparts(filepath);
img_path = fullfile(excel_dir, [excel_base, '_diagnostic.png']);
try
    exportgraphics(fig, img_path, 'Resolution', 200);
    fprintf('Saved constellation diagnostic to:\n  %s\n', img_path);
catch
    % Fallback for older MATLAB versions or if exportgraphics fails.
    try
        saveas(fig, img_path);
        fprintf('Saved constellation diagnostic to:\n  %s\n', img_path);
    catch
        warning('Could not save constellation diagnostic image: %s', img_path);
    end
end

fprintf('Ready to run irradiance model against:\n  %s\n', filepath);

function unique_name = resolve_unique_excel_filename(folder, desired_name)
[~, base_name, ext] = fileparts(desired_name);
if isempty(ext)
    ext = '.xlsx';
end

candidate = [base_name, ext];
serial = 1;
while isfile(fullfile(folder, candidate))
    candidate = sprintf('%s_%03d%s', base_name, serial, ext);
    serial = serial + 1;
end
unique_name = candidate;
end

function update_heliogyro_location_file(location_file_path, excel_file_name)
content = fileread(location_file_path);
pattern = 'excel_file_name\s*=\s*''[^'']*'';';
replacement = sprintf('excel_file_name                    = ''%s'';', excel_file_name);

updated = regexprep(content, pattern, replacement, 'once');
if strcmp(updated, content)
    error('Could not locate excel_file_name assignment in %s', location_file_path);
end

fid = fopen(location_file_path, 'w');
if fid < 0
    error('Failed to open location file for writing: %s', location_file_path);
end
cleanup_obj = onCleanup(@() fclose(fid));
fwrite(fid, updated, 'char');
clear cleanup_obj
fprintf('Updated heliogyro location file to:\n  %s\n', excel_file_name);
end

function render_irradiance_grid_panel(fig, grid_data)
% Replaces the top-left diagnostic subplot with the 11x11 irradiance grid,
% expands it to occupy the full top half, and removes the top-right LD plot.

ax_all = findall(fig, 'Type', 'axes');

% Identify panel 1 and panel 2 by their original titles (from visualize_constellation.m).
ax_grid = [];
ax_ld   = [];
for k = 1:numel(ax_all)
    tObj = get(ax_all(k), 'Title');
    if isempty(tObj) || ~isprop(tObj, 'String')
        continue;
    end
    t = string(tObj.String);
    if isempty(ax_grid) && contains(t, 'Face-on (Y–Z)', 'IgnoreCase', true)
        ax_grid = ax_all(k);
    elseif isempty(ax_ld) && contains(t, 'Solar disk: limb darkening-weighted coverage', 'IgnoreCase', true)
        ax_ld = ax_all(k);
    end
end

% Fallback: if title matching fails, use position heuristics.
if isempty(ax_grid) || isempty(ax_ld)
    axes_pos = arrayfun(@(h) get(h, 'Position'), ax_all, 'UniformOutput', false);
    tops = [];
    for k = 1:numel(ax_all)
        p = axes_pos{k};
        if p(2) > 0.3 % likely top row
            tops(end+1) = k; %#ok<AGROW>
        end
    end
    if numel(tops) >= 2
        ax_candidates = ax_all(tops);
        % left = smaller X position
        pos_candidates = arrayfun(@(h) get(h, 'Position'), ax_candidates, 'UniformOutput', false);
        xs = cellfun(@(p) p(1), pos_candidates);
        [~, order] = sort(xs);
        ax_grid = ax_candidates(order(1));
        ax_ld   = ax_candidates(order(2));
    else
        error('Could not locate expected diagnostic axes in figure.');
    end
end

% Delete the LD coverage panel entirely (including its colorbar, if any).
pos_ld = [];
if ~isempty(ax_ld) && isgraphics(ax_ld)
    pos_ld = get(ax_ld, 'Position');
end
try
    delete(findall(fig, 'Type', 'ColorBar'));
catch
    % no-op (older/newer MATLAB variants)
end
if ~isempty(ax_ld) && isgraphics(ax_ld)
    delete(ax_ld);
end

% Expand the grid axes to span the top half (union of old panel 1 + panel 2 positions).
pos_grid = get(ax_grid, 'Position');
pos_union = pos_grid;
if ~isempty(pos_ld)
    pos_union(3) = (pos_ld(1) + pos_ld(3)) - pos_union(1); % width
end
set(ax_grid, 'Position', pos_union);

% Split the expanded top-half area into:
%   left  -> shading map (in ax_grid)
%   right -> numeric readouts (values at lon=0 and lat=0)
map_w = pos_union(3) * 0.62;
gap_w = pos_union(3) * 0.02;
txt_w = pos_union(3) - map_w - gap_w;

map_pos = [pos_union(1), pos_union(2), map_w, pos_union(4)];
txt_pos = [pos_union(1) + map_w + gap_w, pos_union(2), txt_w, pos_union(4)];

% ----- Panel 1: Earth shading field (normalized lon/lat) -----
cla(ax_grid);
set(ax_grid, 'Position', map_pos);
axes(ax_grid); %#ok<MAXES>

n_rows = size(grid_data, 1);
n_cols = size(grid_data, 2);
latitudes = linspace(-90, 90, n_rows);
longitudes = linspace(-180, 180, n_cols);

lon_norm = longitudes ./ 180;  % [-1, 1]
lat_norm = latitudes  ./ 90;  % [-1, 1]

imagesc(ax_grid, lon_norm, lat_norm, grid_data);
set(ax_grid, 'YDir', 'normal');

colormap(ax_grid, parula(256));
% Keep auto-scaling for contrast.
clim(ax_grid, [min(grid_data(:)), max(grid_data(:))]);

hold(ax_grid, 'on');

% Earth outline (always a true circle in normalized coordinates).
theta = linspace(0, 2*pi, 400);
plot(ax_grid, cos(theta), sin(theta), 'w-', 'LineWidth', 1.6);

% Lightweight lat/lon grid overlays.
grid_lat_lines = [-60 -30 0 30 60];
grid_lon_lines = [0]; % prime meridian
grid_color = [1 1 1] * 0.35;
for la = grid_lat_lines
    plot(ax_grid, [-1 1], [la/90 la/90], '-', 'Color', grid_color, 'LineWidth', 0.6);
end
for lo = grid_lon_lines
    plot(ax_grid, [lo/180 lo/180], [-1 1], '-', 'Color', grid_color, 'LineWidth', 0.6);
end

% Optional continents: keep best-effort and non-fatal.
try
    S = load('coastlines');  % lat/lon (degrees)
    coast_lat = S.lat;
    coast_lon = S.lon;
    coast_lon(coast_lon > 180) = coast_lon(coast_lon > 180) - 360;
    on_disk = (coast_lon/180).^2 + (coast_lat/90).^2 <= 1.02;
    coast_lon_plot = coast_lon;
    coast_lat_plot = coast_lat;
    coast_lon_plot(~on_disk) = NaN;
    coast_lat_plot(~on_disk) = NaN;
    plot(ax_grid, coast_lon_plot/180, coast_lat_plot/90, 'Color', [1 1 1]*0.55, 'LineWidth', 0.6);
catch
    % Skip continents if dataset isn't available.
end

% Formatting.
axis(ax_grid, 'equal');
axis(ax_grid, [-1 1 -1 1]);
set(ax_grid, 'XTick', [-1 -0.5 0 0.5 1], 'YTick', [-1 -0.5 0 0.5 1]);
ax_grid.XTickLabel = arrayfun(@(x) sprintf('%d', x*180), ax_grid.XTick, 'UniformOutput', false);
ax_grid.YTickLabel = arrayfun(@(y) sprintf('%d', y*90),  ax_grid.YTick, 'UniformOutput', false);
grid(ax_grid, 'off');
title(ax_grid, 'Earth shading map', 'FontWeight', 'bold');

% ----- Panel 2: numeric readouts -----
ax_txt = axes('Parent', fig, 'Position', txt_pos);
cla(ax_txt);
axis(ax_txt, 'off');
xlim(ax_txt, [0 1]);
ylim(ax_txt, [0 1]);

mono_font = 'Courier';

% Index for lon=0 column and lat=0 row.
[~, col0_idx] = min(abs(longitudes - 0));
[~, row0_idx] = min(abs(latitudes  - 0));

left_x = 0.05;
right_x = 0.55;

text(ax_txt, left_x, 0.95, 'lon = 0 deg', ...
    'Units', 'normalized', 'FontName', mono_font, 'FontSize', 12, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'Interpreter', 'none');
text(ax_txt, right_x, 0.95, 'lat = 0 deg', ...
    'Units', 'normalized', 'FontName', mono_font, 'FontSize', 12, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'Interpreter', 'none');

% Match analysis_General___0___X area-weighted mean (central longitude column).
central_col = ceil(n_cols / 2);
weights = cosd(latitudes);
weights = weights / sum(weights);
central_strip = grid_data(:, central_col);
mean_irr = sum(central_strip .* weights');

mean_str = sprintf('Area-weighted mean irradiance factor: %.4f', mean_irr);
text(ax_txt, left_x, 0.90, mean_str, ...
    'Units', 'normalized', 'FontName', mono_font, 'FontSize', 12, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'Interpreter', 'none');

top_y = 0.83;
row_step = 0.065;

% Left list: lon = 0 column, labeled by latitude.
for r = 1:n_rows
    y = top_y - (r-1)*row_step;
    lat_val = latitudes(r);
    val = grid_data(r, col0_idx);
    line_str = sprintf('%+6.1f :  %.4f', lat_val, val);
    text(ax_txt, left_x, y, line_str, ...
        'Units', 'normalized', 'FontName', mono_font, 'FontSize', 11, ...
        'HorizontalAlignment', 'left', 'Interpreter', 'none');
end

% Right list: lat = 0 row, labeled by longitude (no repeated "lat/lon =" wording per row).
for c = 1:n_cols
    y = top_y - (c-1)*row_step;
    lon_val = longitudes(c);
    val = grid_data(row0_idx, c);
    line_str = sprintf('%+7.1f :  %.4f', lon_val, val);
    text(ax_txt, right_x, y, line_str, ...
        'Units', 'normalized', 'FontName', mono_font, 'FontSize', 11, ...
        'HorizontalAlignment', 'left', 'Interpreter', 'none');
end
end