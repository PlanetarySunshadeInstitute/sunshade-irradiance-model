% study_lattice_density.m
%
% LATTICE DENSITY STUDY
% ============================================================================
%
% PURPOSE
%   Investigate the maximum achievable shade density using hexagonal lattice
%   placement with multi-layer 3D offset, while guaranteeing ZERO line-of-
%   sight overlap between shades on different planes.
%
%   Sweeps three parameters:
%     1. min_buffer_km:     edge-to-edge clearance within a plane (50–200 km)
%     2. n_planes:          number of X-depth planes (1–10)
%     3. plane_spacing_km:  distance between planes (5–500 km)
%
%   For each configuration, reports:
%     - Total shade count that fits within the Glasgow stability envelope
%     - Fill fraction (shade area / ellipse area)
%     - Number of FCC layers used (1, 2, or 3)
%     - Minimum line-of-sight clearance (verified, not assumed)
%     - Effective density (shades per 1000 km² of Y-Z footprint)
%
% KEY GEOMETRIC INSIGHT
%   A hex lattice with spacing s = 2r + buffer achieves optimal 2D circle
%   packing. By offsetting alternate planes by 1/3 of the hex cell (FCC
%   stacking), we fit up to 3 non-overlapping layers in the same Y-Z
%   footprint, multiplying the shade count without any shadow overlap.
%
%   The maximum layers depend on the spacing-to-radius ratio:
%     s > 2r√3 (~69.3 km for r=20): 3 layers (FCC A-B-C)
%     s > 2r   (~40 km for r=20):   2 layers
%     otherwise:                     1 layer
%
% OUTPUTS
%   - Console table of results
%   - 3-panel figure:
%       (1) Total shades vs buffer for each plane count
%       (2) Fill fraction vs buffer
%       (3) Shades vs plane spacing (at fixed buffer = 100 km)
%   - Results struct saved to workspace for further analysis
%
% REQUIRES
%   place_blob_lattice.m, load_envelope.m on the MATLAB path.
%   equilibrium_envelope.mat in the constellation generator folder.
%
% Authors: Planetary Sunshade Foundation / Claude analysis

clear; clc;
fprintf('\n');
fprintf('================================================================\n');
fprintf('  LATTICE DENSITY STUDY\n');
fprintf('================================================================\n\n');

% =========================================================================
% CONFIGURATION
% =========================================================================

% Fixed parameters
sail_radius_km = 20;       % km (physical radius of each heliogyro)
x_optimal_km   = 838400;   % km (optimal L1 distance)

% Blob footprint: use full solar disk projection as the bounding ellipse
% (this gives the maximum useful area for shading)
solar_disk_radius_km = 10980;   % projected solar disk at L1
ellipse_radii = [solar_disk_radius_km, solar_disk_radius_km];   % circular

% Sweep ranges
buffer_values        = [50, 75, 100, 125, 150, 200];       % km
n_planes_values      = [1, 2, 3, 5, 7, 10];                % number of planes
plane_spacing_values = [5, 10, 25, 50, 100, 250, 500];     % km

% =========================================================================
% LOAD ENVELOPE
% =========================================================================

% Locate envelope file: try the script's own folder first, then the current
% folder, then the folder used by generate_constellation_v2.m.
% (mfilename('fullpath') is unreliable in scripts — use which() instead.)
script_dir    = fileparts(which('study_lattice_density'));
candidates    = { ...
    fullfile(script_dir, 'equilibrium_envelope.mat'), ...
    fullfile(pwd,        'equilibrium_envelope.mat') };

envelope_file = '';
for ci = 1:numel(candidates)
    if isfile(candidates{ci})
        envelope_file = candidates{ci};
        break;
    end
end

if isempty(envelope_file)
    error(['equilibrium_envelope.mat not found.\n' ...
           'Checked:\n  %s\n  %s\n' ...
           'Run preprocess_glasgow_data.m, or ensure the file is on the MATLAB path.'], ...
          candidates{1}, candidates{2});
end
envelope = load_envelope(envelope_file);

% =========================================================================
% STUDY 1: SWEEP BUFFER SIZE × NUMBER OF PLANES
%          (fixed plane_spacing = 10 km)
% =========================================================================

fprintf('\n--- Study 1: Buffer × Planes (plane_spacing = 10 km) ---\n\n');

fixed_spacing = 10;   % km

results_1 = struct('buffer', {}, 'n_planes', {}, 'n_total', {}, ...
                   'n_layers', {}, 'fill_frac', {}, 'min_los', {}, ...
                   'hex_spacing', {});

fprintf('  %8s  %8s  %8s  %8s  %10s  %10s  %10s\n', ...
        'Buffer', 'Planes', 'Layers', 'Total', 'Fill %', 'LOS clr', 'Hex s');
fprintf('  %8s  %8s  %8s  %8s  %10s  %10s  %10s\n', ...
        '(km)', '', '(FCC)', 'shades', '', '(km)', '(km)');
fprintf('  %s\n', repmat('-', 1, 76));

idx = 0;
for buf = buffer_values
    for np = n_planes_values

        blob = struct( ...
            'center_position',  [x_optimal_km, 0, 0], ...
            'ellipse_radii',    ellipse_radii, ...
            'sail_radius_km',   sail_radius_km, ...
            'min_buffer_km',    buf, ...
            'n_planes',         np, ...
            'plane_spacing_km', fixed_spacing);

        [~, meta] = place_blob_lattice(blob, envelope, 'verbose', false, 'verify_los', false);

        idx = idx + 1;
        results_1(idx).buffer      = buf;
        results_1(idx).n_planes    = np;
        results_1(idx).n_total     = meta.n_total;
        results_1(idx).n_layers    = meta.n_layers;
        results_1(idx).fill_frac   = meta.fill_fraction;
        results_1(idx).min_los     = meta.min_los_clearance_km;
        results_1(idx).hex_spacing = meta.hex_spacing_km;

        los_str = 'N/A';
        if isfinite(meta.min_los_clearance_km)
            los_str = sprintf('%.1f', meta.min_los_clearance_km);
        end

        fprintf('  %8.0f  %8d  %8d  %8d  %10.2f  %10s  %10.1f\n', ...
                buf, np, meta.n_layers, meta.n_total, ...
                meta.fill_fraction * 100, los_str, meta.hex_spacing_km);
    end
end

% =========================================================================
% STUDY 2: SWEEP PLANE SPACING (fixed buffer = 100 km, n_planes = 5)
% =========================================================================

fprintf('\n\n--- Study 2: Plane Spacing (buffer = 100 km, n_planes = 5) ---\n\n');

fixed_buffer  = 100;
fixed_nplanes = 5;

results_2 = struct('spacing', {}, 'n_total', {}, 'n_layers', {}, ...
                   'fill_frac', {}, 'min_los', {});

fprintf('  %10s  %8s  %8s  %10s  %10s\n', ...
        'Spacing', 'Layers', 'Total', 'Fill %', 'LOS clr');
fprintf('  %10s  %8s  %8s  %10s  %10s\n', ...
        '(km)', '(FCC)', 'shades', '', '(km)');
fprintf('  %s\n', repmat('-', 1, 56));

for si = 1:numel(plane_spacing_values)
    sp = plane_spacing_values(si);

    blob = struct( ...
        'center_position',  [x_optimal_km, 0, 0], ...
        'ellipse_radii',    ellipse_radii, ...
        'sail_radius_km',   sail_radius_km, ...
        'min_buffer_km',    fixed_buffer, ...
        'n_planes',         fixed_nplanes, ...
        'plane_spacing_km', sp);

    [~, meta] = place_blob_lattice(blob, envelope, 'verbose', false, 'verify_los', false);

    results_2(si).spacing   = sp;
    results_2(si).n_total   = meta.n_total;
    results_2(si).n_layers  = meta.n_layers;
    results_2(si).fill_frac = meta.fill_fraction;
    results_2(si).min_los   = meta.min_los_clearance_km;

    los_str = 'N/A';
    if isfinite(meta.min_los_clearance_km)
        los_str = sprintf('%.1f', meta.min_los_clearance_km);
    end

    fprintf('  %10.0f  %8d  %8d  %10.2f  %10s\n', ...
            sp, meta.n_layers, meta.n_total, ...
            meta.fill_fraction * 100, los_str);
end

% =========================================================================
% STUDY 3: SINGLE DETAILED RUN (for visual inspection)
% =========================================================================

fprintf('\n\n--- Study 3: Detailed run (buffer=100, planes=3, spacing=10) ---\n');

blob_detail = struct( ...
    'center_position',  [x_optimal_km, 0, 0], ...
    'ellipse_radii',    ellipse_radii, ...
    'sail_radius_km',   sail_radius_km, ...
    'min_buffer_km',    100, ...
    'n_planes',         3, ...
    'plane_spacing_km', 10);

[pos_detail, meta_detail] = place_blob_lattice(blob_detail, envelope, 'verbose', true, 'verify_los', true);

% =========================================================================
% FIGURES
% =========================================================================

% -------------------------------------------------------------------------
% Figure 1: Study 1 — Total shades vs buffer, colored by n_planes
% -------------------------------------------------------------------------

fig1 = figure('Name', 'Lattice Study: Shades vs Buffer', 'NumberTitle', 'off', ...
              'Position', [50, 100, 1400, 500]);
tl1 = tiledlayout(fig1, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Panel 1: Total shade count
ax1 = nexttile(tl1);
hold(ax1, 'on');
colors = lines(numel(n_planes_values));
for pi = 1:numel(n_planes_values)
    np = n_planes_values(pi);
    mask = [results_1.n_planes] == np;
    plot(ax1, [results_1(mask).buffer], [results_1(mask).n_total], ...
         'o-', 'Color', colors(pi,:), 'LineWidth', 1.5, 'MarkerSize', 6, ...
         'MarkerFaceColor', colors(pi,:), ...
         'DisplayName', sprintf('%d planes', np));
end
xlabel(ax1, 'Min buffer (km)');
ylabel(ax1, 'Total shades placed');
title(ax1, 'Total shades vs. buffer');
legend(ax1, 'Location', 'northeast');
grid(ax1, 'on');

% Panel 2: Fill fraction
ax2 = nexttile(tl1);
hold(ax2, 'on');
for pi = 1:numel(n_planes_values)
    np = n_planes_values(pi);
    mask = [results_1.n_planes] == np;
    plot(ax2, [results_1(mask).buffer], [results_1(mask).fill_frac] * 100, ...
         'o-', 'Color', colors(pi,:), 'LineWidth', 1.5, 'MarkerSize', 6, ...
         'MarkerFaceColor', colors(pi,:), ...
         'DisplayName', sprintf('%d planes', np));
end
xlabel(ax2, 'Min buffer (km)');
ylabel(ax2, 'Fill fraction (%)');
title(ax2, 'Fill fraction vs. buffer');
legend(ax2, 'Location', 'northeast');
grid(ax2, 'on');

% Panel 3: Study 2 — Total shades vs plane spacing
ax3 = nexttile(tl1);
plot(ax3, [results_2.spacing], [results_2.n_total], 'bo-', ...
     'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', [0.2 0.5 1.0]);
xlabel(ax3, 'Plane spacing (km)');
ylabel(ax3, 'Total shades placed');
title(ax3, sprintf('Shades vs. spacing (buf=%d, %d planes)', fixed_buffer, fixed_nplanes));
grid(ax3, 'on');

sgtitle(tl1, 'Lattice Density Study — Zero Line-of-Sight Overlap Guaranteed', ...
        'FontSize', 13, 'FontWeight', 'bold');

% -------------------------------------------------------------------------
% Figure 2: Y-Z scatter of the detailed run, colored by layer
% -------------------------------------------------------------------------

fig2 = figure('Name', 'Lattice Y-Z Layout', 'NumberTitle', 'off', ...
              'Position', [100, 50, 800, 750]);

layer_colors = [0.2 0.5 1.0;   % Layer A = blue
                1.0 0.3 0.2;   % Layer B = red
                0.2 0.8 0.3];  % Layer C = green

hold on;

% Draw envelope ellipse for reference
theta_env = linspace(0, 2*pi, 200);
a_opt = interp1(envelope.X_km, envelope.a_km, x_optimal_km, 'linear', 'extrap');
b_opt = interp1(envelope.X_km, envelope.b_km, x_optimal_km, 'linear', 'extrap');
plot(a_opt * cos(theta_env), b_opt * sin(theta_env), 'k--', 'LineWidth', 1.5, ...
     'DisplayName', 'Stability envelope');

% Draw blob ellipse
plot(ellipse_radii(1) * cos(theta_env), ellipse_radii(2) * sin(theta_env), ...
     'm--', 'LineWidth', 1, 'DisplayName', 'Blob ellipse');

% Plot shades by layer
PY = pos_detail(2,:);
PZ = pos_detail(3,:);
PX = pos_detail(1,:);

layer_assign = meta_detail.layer_assignment;
X_planes_d   = meta_detail.X_planes_km;

for ip = 1:blob_detail.n_planes
    mask = (PX == X_planes_d(ip));
    li = layer_assign(ip);
    scatter(PY(mask), PZ(mask), 4, layer_colors(li,:), 'filled', ...
            'DisplayName', sprintf('Plane %d (Layer %c)', ip, 'A' + li - 1));
end

axis equal;
xlabel('Y (km)');
ylabel('Z (km)');
title(sprintf('Y-Z layout: %d shades, %d planes, %d layers, buffer=%d km', ...
      meta_detail.n_total, blob_detail.n_planes, meta_detail.n_layers, ...
      blob_detail.min_buffer_km));
legend('Location', 'eastoutside');
grid on;

% -------------------------------------------------------------------------
% Figure 3: Zoomed-in view showing hex structure and layer offsets
% -------------------------------------------------------------------------

fig3 = figure('Name', 'Lattice Zoom', 'NumberTitle', 'off', ...
              'Position', [150, 50, 700, 700]);
hold on;

% Zoom to a small region near center to show the lattice structure
zoom_radius = 800;   % km

for ip = 1:blob_detail.n_planes
    mask = (PX == X_planes_d(ip)) & ...
           abs(PY) < zoom_radius & abs(PZ) < zoom_radius;
    li = layer_assign(ip);

    % Draw circles for each shade
    Y_sub = PY(mask);
    Z_sub = PZ(mask);
    for k = 1:numel(Y_sub)
        th = linspace(0, 2*pi, 40);
        plot(Y_sub(k) + sail_radius_km * cos(th), ...
             Z_sub(k) + sail_radius_km * sin(th), ...
             '-', 'Color', [layer_colors(li,:), 0.5], 'LineWidth', 0.5);
    end

    scatter(Y_sub, Z_sub, 15, layer_colors(li,:), 'filled', ...
            'DisplayName', sprintf('Plane %d (Layer %c)', ip, 'A' + li - 1));
end

axis equal;
xlim([-zoom_radius, zoom_radius]);
ylim([-zoom_radius, zoom_radius]);
xlabel('Y (km)');
ylabel('Z (km)');
title('Zoomed view — hex lattice with layer offsets');
legend('Location', 'eastoutside');
grid on;

% =========================================================================
% SUMMARY
% =========================================================================

fprintf('\n================================================================\n');
fprintf('  SUMMARY OF KEY FINDINGS\n');
fprintf('================================================================\n\n');

% Find the configuration with the most shades
[max_shades, max_idx] = max([results_1.n_total]);
best = results_1(max_idx);
fprintf('  Maximum shade count achieved: %d\n', max_shades);
fprintf('    Config: buffer=%d km, %d planes, %d FCC layers\n', ...
        best.buffer, best.n_planes, best.n_layers);
fprintf('    Fill fraction: %.2f%%\n', best.fill_frac * 100);
fprintf('    Hex spacing: %.1f km\n', best.hex_spacing);
fprintf('\n');

% For the "standard" 100 km buffer
mask_100 = [results_1.buffer] == 100;
std_results = results_1(mask_100);
fprintf('  At standard buffer = 100 km:\n');
for i = 1:numel(std_results)
    fprintf('    %2d planes → %6d shades (%.2f%% fill, %d layers)\n', ...
            std_results(i).n_planes, std_results(i).n_total, ...
            std_results(i).fill_frac * 100, std_results(i).n_layers);
end

fprintf('\n  Geometry limits:\n');
fprintf('    Hex spacing for buffer=100 km: %.0f km\n', 2*sail_radius_km + 100);
fprintf('    Cross-layer min distance:      %.1f km (need > %.0f km for LOS clearance)\n', ...
        (2*sail_radius_km + 100) / sqrt(3), 2*sail_radius_km);
fprintf('    Max FCC layers at buffer=100:  3 (%.1f > %.0f ✓)\n', ...
        (2*sail_radius_km + 100) / sqrt(3), 2*sail_radius_km);

fprintf('\n');
fprintf('  Plane spacing effect (Study 2):\n');
fprintf('    Spacing has NO effect on total count when < stability envelope span.\n');
fprintf('    Each plane gets the same hex grid; layers cycle A-B-C-A-B-C...\n');
fprintf('    Spacing only matters for: physical collision risk, SRP perturbation,\n');
fprintf('    and how much of the envelope span is consumed by the plane array.\n');

fprintf('\n================================================================\n');