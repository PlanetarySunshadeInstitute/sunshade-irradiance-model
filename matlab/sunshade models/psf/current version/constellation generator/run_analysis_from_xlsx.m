% run_analysis_from_xlsx.m
%
% Standalone: point the PSF irradiance analysis at an existing constellation
% .xlsx, refresh the diagnostic figure, and save a .png next to the workbook.
%
% Same post-save sequence as generate_constellation.m (after the user types 'y').
%
% Authors: Planetary Sunshade Foundation

clear; clc;

% -------------------------------------------------------------------------
% USER SETTINGS — edit this block only
% -------------------------------------------------------------------------
FILENAME = 'constellation_polar_11124crafts_2026-03-27_trimodal_v2.xlsx';
% -------------------------------------------------------------------------

filepath = resolve_constellation_xlsx(FILENAME);
[~, fname_stem, fname_ext] = fileparts(filepath);
resolved_name = [fname_stem, fname_ext];

[user_params, xlsx_meta] = parse_constellation_xlsx_headers(filepath);
params = define_constellation_params(user_params);
params.output_filename = resolved_name;

positions = read_constellation_positions_from_xlsx(filepath, params);
if xlsx_meta.refine_polar_radii_from_positions
    params.polar_radius_y_km = max(abs(positions.PY(:)));
    params.polar_radius_z_km = max(abs(positions.PZ(:)));
end
visualize_constellation(positions, params, [], 0);

this_dir = fileparts(mfilename('fullpath'));
location_file = fullfile(this_dir, '..', 'file locations', 'location_Heliogyro_Kinematics_Data___E___Sr.m');
update_heliogyro_location_file(location_file, resolved_name);

results_grid = analysis_General___0___X;
disp('Raw grid (row 1 = first model point, row 11 = last):');
disp(results_grid);
if isempty(results_grid) || ~ismatrix(results_grid)
    error('analysis_General___0___X did not return a 2-D irradiance matrix.');
end
fprintf('\n--- Irradiance factor grid (11x11) ---\n');
disp(round(results_grid, 4));

fig = findobj('Type', 'figure', 'Name', 'Constellation Diagnostic');
if isempty(fig)
    fig = gcf;
else
    fig = fig(1);
end
render_irradiance_grid_panel(fig, results_grid);

img_path = fullfile(fileparts(filepath), [fname_stem, '.png']);
try
    exportgraphics(fig, img_path, 'Resolution', 200);
    fprintf('Saved constellation diagnostic to:\n  %s\n', img_path);
catch
    try
        saveas(fig, img_path);
        fprintf('Saved constellation diagnostic to:\n  %s\n', img_path);
    catch
        warning('Could not save constellation diagnostic image: %s', img_path);
    end
end

fprintf('Ready to run irradiance model against:\n  %s\n', filepath);

% =========================================================================
% Local helpers
% =========================================================================

function filepath = resolve_constellation_xlsx(filename)
if isempty(filename)
    error('FILENAME is empty.');
end
if isfile(filename)
    filepath = filename;
    return
end
folder = location_Excel_Folder___0___St;
candidate = fullfile(folder, filename);
if isfile(candidate)
    filepath = candidate;
    return
end
error('Constellation workbook not found:\n  %s\nnor under Excel folder:\n  %s', ...
      filename, folder);
end

function [user_params, meta] = parse_constellation_xlsx_headers(filepath)
meta = struct('refine_polar_radii_from_positions', false);
metaA = readcell(filepath, 'Sheet', 1, 'Range', 'A1:A9');
lines = cell(9, 1);
for i = 1:9
    if isempty(metaA{i})
        lines{i} = '';
    else
        lines{i} = strtrim(char(string(metaA{i})));
    end
end

tok = regexp(lines{3}, 'Pattern:\s*(.+)', 'tokens', 'once');
if isempty(tok)
    error('Could not parse pattern from row 3 of %s:\n%s', filepath, lines{3});
end
pattern_raw = lower(strtrim(tok{1}));
pattern = normalize_pattern_for_params(pattern_raw);

tok = regexp(lines{4}, 'N craft:\s*(\d+)', 'tokens', 'once');
if isempty(tok)
    error('Could not parse N craft from row 4 of %s:\n%s', filepath, lines{4});
end
N = str2double(tok{1});

tok = regexp(lines{5}, 'N planes:\s*(\d+)\s*\(\s*spacing:\s*([\d.]+)\s*km\s*\)', 'tokens', 'once');
if isempty(tok)
    error('Could not parse N planes / spacing from row 5 of %s:\n%s', filepath, lines{5});
end
n_planes = str2double(tok{1});
plane_spacing_km = str2double(tok{2});

tok = regexp(lines{7}, 'Profile:\s*(\w+)', 'tokens', 'once');
if isempty(tok)
    error('Could not parse profile from row 7 of %s:\n%s', filepath, lines{7});
end
footprint_profile = lower(tok{1});

tok = regexp(lines{8}, 'Sail radius:\s*([\d.]+)\s*km\s*\(\s*min buffer:\s*([\d.]+)\s*km\s*\)', 'tokens', 'once');
if isempty(tok)
    error('Could not parse sail radius / min buffer from row 8 of %s:\n%s', filepath, lines{8});
end
sail_radius_km = str2double(tok{1});
min_buffer_km = str2double(tok{2});

user_params.pattern = pattern;
user_params.N = N;
user_params.n_planes = n_planes;
user_params.plane_spacing_km = plane_spacing_km;
user_params.footprint_profile = footprint_profile;
user_params.sail_radius_km = sail_radius_km;
user_params.min_buffer_km = min_buffer_km;
% Not written to export metadata — match define_constellation_params default.
user_params.footprint_sigma_fraction = 0.4;

if strcmp(pattern, 'uniform')
    tok = regexp(lines{6}, 'Footprint:\s*([\d.]+)\s*km radius\s*\(\s*([\d.]+)%\s*of solar disk\s*\)', 'tokens', 'once');
    if isempty(tok)
        error('Could not parse uniform footprint from row 6 of %s:\n%s', filepath, lines{6});
    end
    user_params.constellation_radius_km = str2double(tok{1});
elseif strcmp(pattern, 'polar')
    fp = lines{6};
    tok_el = regexp(fp, 'Footprint:\s*polar ellipses Y\s*([\d.]+)\s*km,\s*Z\s*([\d.]+)\s*km', 'tokens', 'once');
    tok_tri = regexp(fp, '(?i)low\s*Z\s*(-?\d+)\s*,\s*mid\s*Z\s*(-?\d+)\s*,\s*high\s*Z\s*(-?\d+)', 'tokens', 'once');
    if ~isempty(tok_el)
        user_params.polar_radius_y_km = str2double(tok_el{1});
        user_params.polar_radius_z_km = str2double(tok_el{2});
        polar = read_polar_block_after_data(filepath, N);
        user_params.polar_center_z_km = polar.center_z_km;
        user_params.hemisphere = polar.hemisphere;
        user_params.polar_fraction = polar.fraction;
    elseif ~isempty(tok_tri)
        zl = str2double(tok_tri{1});
        zm = str2double(tok_tri{2});
        zh = str2double(tok_tri{3});
        user_params.polar_center_z_km = mean([zl, zm, zh]);
        user_params.hemisphere = 'both';
        user_params.polar_fraction = 1;
        user_params.polar_radius_y_km = 2000;
        user_params.polar_radius_z_km = 2000;
        meta.refine_polar_radii_from_positions = true;
    else
        user_params.polar_center_z_km = 6200;
        user_params.hemisphere = 'both';
        user_params.polar_fraction = 1;
        user_params.polar_radius_y_km = 2000;
        user_params.polar_radius_z_km = 2000;
        meta.refine_polar_radii_from_positions = true;
        warning('run_analysis_from_xlsx:FootprintText', ...
            ['Footprint row is not standard polar or trimodal; ', ...
             'using placeholder Z center and ellipse axes from sail positions.\n%s'], fp);
        polar = read_polar_block_after_data(filepath, N);
        if ~isnan(polar.center_z_km)
            user_params.polar_center_z_km = polar.center_z_km;
            user_params.hemisphere = polar.hemisphere;
            user_params.polar_fraction = polar.fraction;
        end
    end
else
    error('Unsupported pattern "%s" (normalized from "%s") in %s', ...
          pattern, pattern_raw, filepath);
end
end

function pattern = normalize_pattern_for_params(pattern_raw)
if isempty(pattern_raw)
    error('Empty pattern string.');
end
if strcmp(pattern_raw, 'uniform') || strcmp(pattern_raw, 'polar')
    pattern = pattern_raw;
    return
end
if startsWith(pattern_raw, 'polar_blob') || contains(pattern_raw, 'blob_opt') || ...
        contains(pattern_raw, 'polar_blob')
    pattern = 'polar';
    return
end
error(['Unsupported constellation pattern "%s". ', ...
       'Expected uniform, polar, or polar_blob* / optimizer id.'], pattern_raw);
end

function polar = read_polar_block_after_data(filepath, N)
polar = struct('center_z_km', NaN, 'hemisphere', 'both', 'fraction', NaN);
frac_found = false;
start_row = N + 12;
blk = readcell(filepath, 'Sheet', 1, 'Range', sprintf('A%d:A%d', start_row, start_row + 8));
for i = 1:numel(blk)
    if isempty(blk{i})
        continue
    end
    s = strtrim(char(string(blk{i})));
    tok = regexp(s, 'Polar center Z:\s*([\d.]+)\s*km\s*\(([^)]+)\)', 'tokens', 'once');
    if ~isempty(tok)
        polar.center_z_km = str2double(tok{1});
        polar.hemisphere = lower(strtrim(tok{2}));
    end
    tok = regexp(s, 'Polar fraction:\s*([\d.]+)%', 'tokens', 'once');
    if ~isempty(tok)
        polar.fraction = str2double(tok{1}) / 100;
        frac_found = true;
    end
end
if isnan(polar.center_z_km)
    warning(['Polar metadata block (rows after data) missing or unparseable — ', ...
             'using defaults for polar_center_z_km / hemisphere / polar_fraction.']);
    polar.center_z_km = 6200;
    polar.hemisphere = 'both';
    polar.fraction = 0.33;
elseif ~frac_found
    polar.fraction = 1;
end
end

function positions = read_constellation_positions_from_xlsx(filepath, params)
mat = readmatrix(filepath, 'Sheet', 1, 'Range', 'A11:C20000');
if size(mat, 2) ~= 3
    error('Expected three columns (PX,PY,PZ) from A11 in %s', filepath);
end
good = ~any(isnan(mat), 2);
mat = mat(good, :);
n = size(mat, 1);
if n ~= params.N
    error('Row count after NaN removal (%d) does not match header N craft (%d) in %s.', ...
          n, params.N, filepath);
end

positions.PX = mat(:, 1);
positions.PY = mat(:, 2);
positions.PZ = mat(:, 3);
positions.NX = ones(n, 1);
positions.NY = zeros(n, 1);
positions.NZ = zeros(n, 1);
positions.pattern = params.pattern;

Xp = params.X_planes_km(:)';
plane_id = zeros(n, 1);
for i = 1:n
    [~, idx] = min(abs(positions.PX(i) - Xp));
    plane_id(i) = idx;
end
positions.plane_id = plane_id;
end

function update_heliogyro_location_file(location_file_path, excel_file_name)
content = fileread(location_file_path);
if isempty(content)
    error('Empty or unreadable file: %s', location_file_path);
end
if contains(content, sprintf('\r\n'))
    brk = sprintf('\r\n');
else
    brk = sprintf('\n');
end
lines = regexp(content, '\r?\n', 'split');
replacement_body = sprintf('excel_file_name                    = ''%s'';', excel_file_name);
hit = false;
for k = 1:numel(lines)
    if isempty(regexp(lines{k}, '^\s*excel_file_name\s*=', 'once'))
        continue
    end
    m = regexp(lines{k}, '^(\s*)', 'tokens', 'once');
    indent = '';
    if ~isempty(m)
        indent = m{1};
    end
    lines{k} = [indent, replacement_body]; %#ok<AGROW>
    hit = true;
    break
end
if ~hit
    error('Could not locate excel_file_name assignment in %s', location_file_path);
end
updated = strjoin(lines, brk);

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
ax_all = findall(fig, 'Type', 'axes');

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

if isempty(ax_grid) || isempty(ax_ld)
    axes_pos = arrayfun(@(h) get(h, 'Position'), ax_all, 'UniformOutput', false);
    tops = [];
    for k = 1:numel(ax_all)
        p = axes_pos{k};
        if p(2) > 0.3
            tops(end+1) = k; %#ok<AGROW>
        end
    end
    if numel(tops) >= 2
        ax_candidates = ax_all(tops);
        pos_candidates = arrayfun(@(h) get(h, 'Position'), ax_candidates, 'UniformOutput', false);
        xs = cellfun(@(p) p(1), pos_candidates);
        [~, order] = sort(xs);
        ax_grid = ax_candidates(order(1));
        ax_ld   = ax_candidates(order(2));
    else
        error('Could not locate expected diagnostic axes in figure.');
    end
end

pos_ld = [];
if ~isempty(ax_ld) && isgraphics(ax_ld)
    pos_ld = get(ax_ld, 'Position');
end
try
    delete(findall(fig, 'Type', 'ColorBar'));
catch
end
if ~isempty(ax_ld) && isgraphics(ax_ld)
    delete(ax_ld);
end

pos_grid = get(ax_grid, 'Position');
pos_union = pos_grid;
if ~isempty(pos_ld)
    pos_union(3) = (pos_ld(1) + pos_ld(3)) - pos_union(1);
end
set(ax_grid, 'Position', pos_union);

map_w = pos_union(3) * 0.62;
gap_w = pos_union(3) * 0.02;
txt_w = pos_union(3) - map_w - gap_w;

map_pos = [pos_union(1), pos_union(2), map_w, pos_union(4)];
txt_pos = [pos_union(1) + map_w + gap_w, pos_union(2), txt_w, pos_union(4)];

cla(ax_grid);
set(ax_grid, 'Position', map_pos);
axes(ax_grid); %#ok<MAXES>

n_rows = size(grid_data, 1);
n_cols = size(grid_data, 2);
latitudes = linspace(-90, 90, n_rows);
longitudes = linspace(-180, 180, n_cols);

lon_norm = longitudes ./ 180;
lat_norm = latitudes  ./ 90;

imagesc(ax_grid, lon_norm, lat_norm, grid_data);
set(ax_grid, 'YDir', 'normal');

colormap(ax_grid, parula(256));
clim(ax_grid, [min(grid_data(:)), max(grid_data(:))]);

hold(ax_grid, 'on');

theta = linspace(0, 2*pi, 400);
plot(ax_grid, cos(theta), sin(theta), 'w-', 'LineWidth', 1.6);

grid_lat_lines = [-60 -30 0 30 60];
grid_lon_lines = [0];
grid_color = [1 1 1] * 0.35;
for la = grid_lat_lines
    plot(ax_grid, [-1 1], [la/90 la/90], '-', 'Color', grid_color, 'LineWidth', 0.6);
end
for lo = grid_lon_lines
    plot(ax_grid, [lo/180 lo/180], [-1 1], '-', 'Color', grid_color, 'LineWidth', 0.6);
end

try
    S = load('coastlines');
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
end

axis(ax_grid, 'equal');
axis(ax_grid, [-1 1 -1 1]);
set(ax_grid, 'XTick', [-1 -0.5 0 0.5 1], 'YTick', [-1 -0.5 0 0.5 1]);
ax_grid.XTickLabel = arrayfun(@(x) sprintf('%d', x*180), ax_grid.XTick, 'UniformOutput', false);
ax_grid.YTickLabel = arrayfun(@(y) sprintf('%d', y*90),  ax_grid.YTick, 'UniformOutput', false);
grid(ax_grid, 'off');
title(ax_grid, 'Earth shading map', 'FontWeight', 'bold');

ax_txt = axes('Parent', fig, 'Position', txt_pos);
cla(ax_txt);
axis(ax_txt, 'off');
xlim(ax_txt, [0 1]);
ylim(ax_txt, [0 1]);

mono_font = 'Courier';

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

for r = 1:n_rows
    y = top_y - (r-1)*row_step;
    lat_val = latitudes(r);
    val = grid_data(r, col0_idx);
    line_str = sprintf('%+6.1f :  %.4f', lat_val, val);
    text(ax_txt, left_x, y, line_str, ...
        'Units', 'normalized', 'FontName', mono_font, 'FontSize', 11, ...
        'HorizontalAlignment', 'left', 'Interpreter', 'none');
end

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
