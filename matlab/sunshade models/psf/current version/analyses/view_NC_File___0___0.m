
%% view_NC_File___0___0.m
%  Seasonal irradiance-factor viewer — 2x2 imagesc diagnostic.
%
%  Reads a preconfigured-mode NC file and shows the sun-facing hemisphere
%  irradiance factor at four key dates (equinoxes and solstices) in a
%  2x2 grid.  A statistics annotation lists file-level metadata.
%
%  Coordinate convention
%  ---------------------
%    Rows  (Y) : phi   — degrees from sub-stellar point, vertical   (−90° to +90°)
%    Columns (X) : theta — degrees from sub-stellar point, horizontal (−90° to +90°)
%
%  The NC file stores data in CESM lon/lat order (288 lon × 192 lat × 365 days).
%  This script reconstructs the 144-column sun-facing hemisphere by reversing
%  the lon-mapping applied in analysis_General___0___X:
%    CESM lon cols   1– 72  (  0°– 88.75°) → disc theta +0.625° to +89.375°  (east)
%    CESM lon cols 217–288  (270°–358.75°) → disc theta −89.375° to −0.625°  (west)
%
%  Usage
%  -----
%    1. Set nc_file_path to the NC file you want to inspect.
%    2. Optionally adjust time_indices to select different days of year.
%    3. Run the script.  A figure opens and is saved as a PNG alongside the NC file.
%
%  Output
%  ------
%    Interactive figure + PNG saved in the same folder as the NC file.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    USER SETTINGS                                                                                                                                            %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   NC file to inspect.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


nc_file_path  = '/Users/morgangoodwin/Desktop/PSI/Matlab/numerical control/exports/L1-RX-2026-04-17-001.nc';       % ← set full path to your NC file here, e.g.:
                          %   '/path/to/exports/L1-RX-2026-04-17-001.nc'


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Day-of-year indices for the four panels.
%   Defaults: March equinox (~80), June solstice (~172),
%             September equinox (~266), December solstice (~355).
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


time_indices  = [80, 172, 266, 355];
panel_labels  = {'March equinox', 'June solstice', 'September equinox', 'December solstice'};


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Set to true to save a PNG alongside the NC file; false to leave the figure open only.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


save_png      = true;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Gamma exponent for colormap stretch.
%   gamma < 1 allocates more colors to the shaded (low IF, blue) end of the
%   range, making subtle gradients in the shaded region more visible.
%   gamma = 1.0 is a linear colormap (no stretch).
%   Suggested range: 0.2 – 0.5 for typical polar constellation files.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


gamma         = 0.35;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    VALIDATE INPUT                                                                                                                                           %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Check the NC file path was set.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


if isempty(nc_file_path)
    error('view_NC_File:MissingPath', ...
        'Set nc_file_path to the full path of your NC file before running.');
end

if ~isfile(nc_file_path)
    error('view_NC_File:FileNotFound', ...
        'NC file not found:\n  %s', nc_file_path);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    READ NC FILE                                                                                                                                             %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Read coordinate vectors and the irradiance-factor array.
%   df is stored as [lon x lat x time] = [288 x 192 x 365].
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


fprintf('Reading NC file...\n  %s\n', nc_file_path);

phi_deg      = ncread(nc_file_path, 'lat');    % 192 x 1, phi values (sub-stellar degrees, vertical)
df_full      = ncread(nc_file_path, 'df');     % 288 x 192 x 365


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Read date strings for the selected time steps (5-char strings: 'MM-DD').
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


date_chars   = ncread(nc_file_path, 'date');   % string_len x n_days


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Validate time indices against file size.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


n_times = size(df_full, 3);
if any(time_indices > n_times) || any(time_indices < 1)
    error('view_NC_File:BadTimeIndex', ...
        'time_indices must be in [1, %d] for this file.', n_times);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    RECONSTRUCT SUN-FACING HEMISPHERE                                                                                                                       %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Re-order the 288 CESM lon columns back into the 144-column sub-stellar disc.
%   Result: disc_idx maps CESM lon → disc theta order, west to east.
%
%   CESM lon cols 217:288 → disc theta −89.375° to −0.625°  (western dayside)
%   CESM lon cols   1: 72 → disc theta  +0.625° to +89.375° (eastern dayside)
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


disc_idx     = [217:288, 1:72];                     % 144 CESM lon indices
n_theta      = numel(disc_idx);                     % 144
n_phi        = numel(phi_deg);                      % 192
theta_deg    = linspace(-89.375, 89.375, n_theta);  % disc theta axis


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Extract the four panels as [phi x theta] matrices (rows = phi, cols = theta).
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


panels = zeros(n_phi, n_theta, 4);

for k = 1 : 4
    t         = time_indices(k);
    slice_lon = squeeze(df_full(disc_idx, :, t));   % 144 x 192
    panels(:, :, k) = slice_lon';                   % 192 x 144  (phi x theta)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    COLORMAP AND SCALE                                                                                                                                       %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Scale the colormap to the actual range across all four panels.
%   Blue = most shading (lowest irradiance factor); white = no shading (1.0).
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


v_min = min(panels(:));
v_max = max(panels(:));

if abs(v_max - v_min) < 1e-8
    v_max = v_min + 1e-4;   % guard against flat data
end

%   frac_lin  : linear ramp over [0, 1] — maps colormap index to data position.
%   frac_g    : gamma-stretched ramp — compresses the white (high IF) end and
%               expands the blue (low IF, shaded) end so more distinct colors
%               fall where the shading gradient actually is.
%   The clim and all data values are unchanged; only the color allocation shifts.


n_colors  = 256;
frac_lin  = linspace(0, 1, n_colors)';
frac_g    = frac_lin .^ gamma;

r_ch      = interp1([0, 1], [0.10, 1.00], frac_g);
g_ch      = interp1([0, 1], [0.30, 1.00], frac_g);
b_ch      = interp1([0, 1], [0.80, 1.00], frac_g);
cmap      = [r_ch, g_ch, b_ch];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    STATISTICS                                                                                                                                               %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Area-weighted mean irradiance factor over the sun-facing hemisphere,
%   computed for each of the four panels.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


weights      = cosd(phi_deg(:));
weights      = weights / sum(weights);

panel_means  = zeros(1, 4);
for k = 1 : 4
    col_means       = weights' * panels(:, :, k);   % 1 x n_theta
    panel_means(k)  = mean(col_means);
end

global_mean  = mean(panel_means);
frac_shaded  = mean(panels(:) < (v_max - 1e-6)) * 100;   % fraction of cells below peak


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Build the statistics annotation string.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


[~, nc_stem, ~]  = fileparts(nc_file_path);

stats_str = sprintf([ ...
    'File:  %s\n',                           ...
    'Grid:  %d phi  x  %d theta\n',          ...
    'Time:  %d-day file,  4 samples\n',      ...
    '\n',                                    ...
    'Dimming factor\n',                      ...
    '  Min          %.8f\n',                 ...
    '  Max          %.8f\n',                 ...
    '  Range        %.2e\n',                 ...
    '  Global mean  %.8f\n',                 ...
    '  Cells < max  %.2f%%\n',               ...
    '\n',                                    ...
    'Panel means\n',                         ...
    '  Mar equinox  %.8f\n',                 ...
    '  Jun solstice %.8f\n',                 ...
    '  Sep equinox  %.8f\n',                 ...
    '  Dec solstice %.8f\n',                 ...
    ], nc_stem, n_phi, n_theta, n_times, ...
    v_min, v_max, v_max - v_min, global_mean, frac_shaded, ...
    panel_means(1), panel_means(2), panel_means(3), panel_means(4));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    FIGURE                                                                                                                                                   %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Create figure: 2x2 imagesc grid + statistics panel on the right.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


fig = figure( ...
    'Name',        'NC File Viewer — Dimming Factor', ...
    'NumberTitle', 'off', ...
    'Color',       'white', ...
    'Position',    [80, 80, 1300, 820] );

tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
tl.Title.String   = sprintf('Dimming factor — sun-facing hemisphere   |   %s', nc_stem);
tl.Title.FontSize = 12;
tl.Title.FontWeight = 'bold';


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Panel positions: tiles 1-4 are the 2x2 imagesc grid; tile 5 (col 3, rows 1-2) is statistics.
%   We get the 2x2 in columns 1-2 by spanning two rows in a 2x3 layout.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


panel_positions = [1, 2, 4, 5];   % tiles in a 2x3 layout: top-left, top-mid, bot-left, bot-mid


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Draw the four panels.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


ax_panels = gobjects(4, 1);

for k = 1 : 4


    ax = nexttile(tl, panel_positions(k));
    ax_panels(k) = ax;

    imagesc(ax, theta_deg, phi_deg, panels(:, :, k));
    set(ax, 'YDir', 'normal');

    colormap(ax, cmap);
    clim(ax, [v_min, v_max]);

    %   Read the date string from the NC file for this time index.
    try
        date_str = strtrim(char(date_chars(:, time_indices(k))'));
    catch
        date_str = sprintf('day %d', time_indices(k));
    end

    title(ax, sprintf('%s  (day %d,  %s)', panel_labels{k}, time_indices(k), date_str), ...
        'FontSize', 10, 'FontWeight', 'bold');

    xlabel(ax, 'Theta — degrees from sub-stellar point [deg]', 'FontSize', 9);
    ylabel(ax, 'Phi — degrees from sub-stellar point [deg]',   'FontSize', 9);

    xlim(ax, [theta_deg(1), theta_deg(end)]);
    ylim(ax, [phi_deg(1),   phi_deg(end)  ]);

    xticks(ax,      [-90, -60, -30, 0, 30, 60, 90]);
    xticklabels(ax, {'-90°','-60°','-30°','0°','30°','60°','90°'});
    yticks(ax,      [-90, -60, -30, 0, 30, 60, 90]);
    yticklabels(ax, {'-90°','-60°','-30°','0°','30°','60°','90°'});

    set(ax, 'TickDir', 'out', 'Box', 'on', 'FontSize', 9, 'Layer', 'top');

    %   Grid lines at ±30° and ±60° in both theta (horizontal) and phi (vertical).
    %   'Layer','top' on the axes ensures these render above the imagesc data.
    hold(ax, 'on');
    grid_vals = [-60, -30, 0, 30, 60];
    for gv = grid_vals
        plot(ax, [theta_deg(1), theta_deg(end)], [gv, gv], ...
             'Color', [0.5 0.5 0.5], 'LineWidth', 0.4, 'LineStyle', '--');
        plot(ax, [gv, gv], [phi_deg(1), phi_deg(end)], ...
             'Color', [0.5 0.5 0.5], 'LineWidth', 0.4, 'LineStyle', '--');
    end
    hold(ax, 'off');

end


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Shared colorbar spanning column 1-2.  Placed explicitly inside the last panel.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


cb              = colorbar(ax_panels(4), 'southoutside');
cb.Label.String = sprintf('Dimming factor   [%.6f  –  %.6f]   (colormap gamma = %.2f)', v_min, v_max, gamma);
cb.Label.FontSize = 9;
cb.FontSize       = 8;
cb.TickLabelsMode = 'manual';
cb.TickLabels     = arrayfun(@(v) sprintf('%.6f', v), cb.Ticks, 'UniformOutput', false);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Statistics panel in tile 3 and 6 (right column, both rows).
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


ax_stats = nexttile(tl, 3, [2, 1]);
axis(ax_stats, 'off');
title(ax_stats, 'Statistics', 'FontSize', 11, 'FontWeight', 'bold', 'Visible', 'on');

text(ax_stats, 0.04, 0.92, stats_str, ...
    'Units',               'normalized', ...
    'VerticalAlignment',   'top', ...
    'HorizontalAlignment', 'left', ...
    'FontSize',            9, ...
    'FontName',            'Courier New', ...
    'Interpreter',         'none');


drawnow;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    SAVE                                                                                                                                                     %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Save the figure as a PNG alongside the NC file.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


if save_png

    [nc_folder, nc_stem_save, ~] = fileparts(nc_file_path);
    png_path = fullfile(char(nc_folder), sprintf('%s_seasonal_viewer.png', char(nc_stem_save)));

    try
        exportgraphics(fig, png_path, 'Resolution', 300);
        fprintf('Saved viewer PNG:\n  %s\n', png_path);
    catch
        try
            saveas(fig, png_path, 'png');
            fprintf('Saved viewer PNG (saveas):\n  %s\n', png_path);
        catch err
            warning('Could not save PNG: %s', err.message);
        end
    end

    %   Also save to writeup/figures/ for the methods paper.
    %   Path is resolved relative to the location of this script file.
    repo_root      = fileparts(fileparts(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))))));
    writeup_path   = fullfile(repo_root, 'writeup', 'figures', 'fig_seasonal_disc.png');
    writeup_dir    = fileparts(writeup_path);
    if ~isfolder(writeup_dir), mkdir(writeup_dir); end
    try
        exportgraphics(fig, writeup_path, 'Resolution', 300);
        fprintf('Saved writeup figure:\n  %s\n', writeup_path);
    catch err
        warning('Could not save writeup figure: %s', err.message);
    end

end


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
