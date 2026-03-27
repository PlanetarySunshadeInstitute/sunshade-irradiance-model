%% plot_shading_map_earth.m
% Standalone diagnostic: 2D Earth shading map from an 11x11 irradiance-factor grid.
%
% Expected input:
%   grid_data : 11x11 matrix of irradiance factors
%               rows    -> latitude (deg)  from -90 to 90
%               columns -> longitude (deg) from -180 to 180
%
% What it produces:
%   Figure with two panels:
%     [1] Shading field (color map) + Earth disk boundary + simple lat/lon grid
%     [2] Readout lists for lon=0 column and lat=0 row (4 decimals, monospaced)
%
% Save:
%   Uses exportgraphics to write a PNG.
%
% Usage:
%   1) In MATLAB, create `grid_data` (11x11) in the base workspace.
%   2) Run this script.

% Do not clear the workspace: the 11x11 grid is expected to be produced
% by the upstream analysis script and may have a different variable name.
clc;

% -------------------------------
% USER SETTINGS
% -------------------------------
output_png = ''; % leave '' to auto-place next to the running script

% If you want to override colormap scaling:
use_fixed_clim = false;
fixed_clim_min = 0.9900;
fixed_clim_max = 1.0000;

% -------------------------------
% VALIDATE INPUT
% -------------------------------
if ~exist('grid_data', 'var') || isempty(grid_data)
    % Be forgiving about the variable name: look for a few common candidates.
    candidate_names = {'grid_data', 'results_grid', 'irradiance_grid', 'irradiance_factor_grid'};
    found = false;
    for i = 1:numel(candidate_names)
        nm = candidate_names{i};
        if evalin('base', sprintf('exist(''%s'',''var'')', nm))
            tmp = evalin('base', nm);
            if isnumeric(tmp) && ismatrix(tmp) && all(size(tmp) == [11 11])
                grid_data = tmp; %#ok<NASGU>
                found = true;
                break;
            end
        end
    end

    % Final fallback: pick the first numeric 11x11 variable in the base workspace.
    if ~found
        ws = evalin('base', 'whos');
        for k = 1:numel(ws)
            if isnumeric(evalin('base', ws(k).name)) && ismatrix(evalin('base', ws(k).name)) && all(ws(k).size == [11 11])
                grid_data = evalin('base', ws(k).name);
                fprintf('[plot_shading_map_earth] Using base workspace variable: %s\n', ws(k).name);
                found = true;
                break;
            end
        end
    end

    if ~found
        error('plot_shading_map_earth:MissingInput', ...
            'Provide an 11x11 irradiance-factor matrix as `grid_data` (or another 11x11 numeric matrix in the base workspace) before running.');
    end
end

if ~ismatrix(grid_data) || ~all(size(grid_data) == [11 11])
    error('plot_shading_map_earth:BadSize', ...
        '`grid_data` must be an 11x11 matrix. Current size: %dx%d', size(grid_data, 1), size(grid_data, 2));
end

% -------------------------------
% COORDINATES (match 11x11 grid spec)
% -------------------------------
n_rows = size(grid_data, 1);
n_cols = size(grid_data, 2);

latitudes_deg  = linspace(-90, 90, n_rows);
longitudes_deg = linspace(-180, 180, n_cols);

% -------------------------------
% FIGURE SETUP
% -------------------------------
fig = figure('Color', 'w', 'Name', 'Earth Shading Map', ...
              'NumberTitle', 'off', 'Position', [120 80 1100 650]);
clf(fig);

% Panel 1: shading field
ax_map = subplot(1, 2, 1, 'Parent', fig);
cla(ax_map);
axes(ax_map);

% Plot in normalized lon/lat so the Earth outline is always a circle,
% independent of figure aspect ratio.
lon_norm = longitudes_deg ./ 180;  % in [-1, 1]
lat_norm = latitudes_deg  ./ 90;   % in [-1, 1]

% Color field: low irradiance (more shading) = blue; high irradiance (less shading) = warm.
imagesc(ax_map, lon_norm, lat_norm, grid_data);
set(ax_map, 'YDir', 'normal');

% Colormap: default MATLAB 'parula' is blue->yellow (warm). Matches the "blue more shading" requirement.
colormap(ax_map, parula(256));

if use_fixed_clim
    clim(ax_map, [fixed_clim_min, fixed_clim_max]);
else
    % Auto-scale for visible contrast while keeping values comparable.
    clim(ax_map, [min(grid_data(:)), max(grid_data(:))]);
end

cb = colorbar(ax_map);
cb.Label.String = 'Irradiance factor';
cb.FontSize = 10;

hold(ax_map, 'on');

% -------------------------------
% Earth disk boundary overlay
% -------------------------------
% Use a simple circle in lon/lat axes coordinates:
%   (lon/180)^2 + (lat/90)^2 = 1
theta = linspace(0, 2*pi, 400);
earth_lon_norm = cos(theta);
earth_lat_norm = sin(theta);
plot(ax_map, earth_lon_norm, earth_lat_norm, 'w-', 'LineWidth', 1.6);

% -------------------------------
% Latitude/longitude grid overlay
% -------------------------------
grid_lat_lines = [-60 -30 0 30 60];
grid_lon_lines = [0]; % prime meridian

grid_color = [1 1 1] * 0.35; % light/low-contrast
for la = grid_lat_lines
    la_norm = la ./ 90;
    plot(ax_map, [-1 1], [la_norm la_norm], '-', 'Color', grid_color, 'LineWidth', 0.6);
end
for lo = grid_lon_lines
    lo_norm = lo ./ 180;
    plot(ax_map, [lo_norm lo_norm], [-1 1], '-', 'Color', grid_color, 'LineWidth', 0.6);
end

% -------------------------------
% Simple continent outline (optional)
% -------------------------------
% If Mapping Toolbox is available, you could use geoshow/worldmap with a low-res land layer.
% To keep this script self-contained, we try MATLAB's built-in `coastlines` dataset.
try
    S = load('coastlines');  % provides variables `lat` and `lon` (degrees)
    coast_lat = S.lat;
    coast_lon = S.lon;

    % Ensure longitude is in [-180, 180]
    coast_lon(coast_lon > 180) = coast_lon(coast_lon > 180) - 360;

    % Mask points outside the disk to avoid distracting lines.
    on_disk = (coast_lon/180).^2 + (coast_lat/90).^2 <= 1.02;
    coast_lon_plot = coast_lon;
    coast_lat_plot = coast_lat;
    coast_lon_plot(~on_disk) = NaN;
    coast_lat_plot(~on_disk) = NaN;

    % Plot in normalized coordinates to match the circle overlay.
    plot(ax_map, coast_lon_plot ./ 180, coast_lat_plot ./ 90, ...
        'Color', [1 1 1] * 0.55, 'LineWidth', 0.6);
catch
    % Skip continents if the dataset isn't available.
end

% -------------------------------
% Formatting
% -------------------------------
axis(ax_map, 'tight');
xlabel(ax_map, 'Longitude (deg)', 'FontSize', 12);
ylabel(ax_map, 'Latitude (deg)', 'FontSize', 12);
title(ax_map, 'Earth shading map', 'FontWeight', 'bold', 'FontSize', 14);
% Force square scaling in normalized space so the Earth outline renders as a true circle.
axis(ax_map, [-1 1 -1 1]);
axis(ax_map, 'equal');

% Keep ticks labeled in degrees for readability.
xt = linspace(-1, 1, 5);
yt = linspace(-1, 1, 5);
ax_map.XTick = xt;
ax_map.YTick = yt;
ax_map.XTickLabel = arrayfun(@(x) sprintf('%d', round(x * 180)), xt, 'UniformOutput', false);
ax_map.YTickLabel = arrayfun(@(y) sprintf('%d', round(y * 90)), yt, 'UniformOutput', false);

% -------------------------------
% Panel 2: value readout lists
% -------------------------------
ax_txt = subplot(1, 2, 2, 'Parent', fig);
cla(ax_txt);
axis(ax_txt, 'off');
xlim(ax_txt, [0 1]);
ylim(ax_txt, [0 1]);

mono_font = 'Courier';
readable_font = 12;

% Index for lon=0 column and lat=0 row.
[~, col0_idx] = min(abs(longitudes_deg - 0));
[~, row0_idx] = min(abs(latitudes_deg - 0));

left_x = 0.05;
right_x = 0.55;

text(ax_txt, left_x, 0.95, 'longitude', ...
    'Units', 'normalized', 'FontName', mono_font, 'FontSize', readable_font, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'Interpreter', 'none');
text(ax_txt, right_x, 0.95, 'latitude', ...
    'Units', 'normalized', 'FontName', mono_font, 'FontSize', readable_font, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'Interpreter', 'none');

top_y = 0.88;
row_step = 0.065;

% Left list: lon=0 column, labeled by latitude
for r = 1:n_rows
    y = top_y - (r-1)*row_step;
    lat_val = latitudes_deg(r);
    val = grid_data(r, col0_idx);
    line_str = sprintf('%+6.1f :  %.4f', lat_val, val);
    text(ax_txt, left_x, y, line_str, ...
        'Units', 'normalized', ...
        'FontName', mono_font, 'FontSize', 11, ...
        'HorizontalAlignment', 'left', 'Interpreter', 'none');
end

% Right list: lat=0 row, labeled by longitude
for c = 1:n_cols
    y = top_y - (c-1)*row_step;
    lon_val = longitudes_deg(c);
    val = grid_data(row0_idx, c);
    line_str = sprintf('%+7.1f :  %.4f', lon_val, val);
    text(ax_txt, right_x, y, line_str, ...
        'Units', 'normalized', ...
        'FontName', mono_font, 'FontSize', 11, ...
        'HorizontalAlignment', 'left', 'Interpreter', 'none');
end

% -------------------------------
% Save
% -------------------------------
if isempty(output_png)
    % Auto-place next to this script (best-effort).
    script_path = mfilename('fullpath');
    [script_dir, ~, ~] = fileparts(script_path);
    output_png = fullfile(script_dir, 'earth_shading_map_diagnostic.png');
end

try
    exportgraphics(fig, output_png, 'Resolution', 300);
    fprintf('Saved Earth shading map to:\n  %s\n', output_png);
catch ME
    warning('Could not save via exportgraphics: %s', ME.message);
end

