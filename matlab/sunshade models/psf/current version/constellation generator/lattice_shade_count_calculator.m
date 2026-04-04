% lattice_shade_count_calculator.m
%
% INVERSE LATTICE CALCULATOR
% ============================================================================
%
% PURPOSE
%   Work backward from a desired shade count to find the blob ellipse radii
%   that produce it. This is the inverse of what place_blob_lattice does:
%   instead of (radii → count), we compute (target count → radii).
%
% BACKGROUND: LATTICE GEOMETRY
%   The hexagonal lattice places shades on a hex grid with center-to-center
%   spacing s = 2*sail_radius + min_buffer. The area of one hex cell is:
%
%       A_cell = s² × √3/2   [km²]
%
%   The number of shades in an ellipse of semi-axes (Ry, Rz) on a single
%   hex layer is approximately:
%
%       N_per_plane ≈ (π × Ry × Rz) / A_cell
%
%   With n_planes planes (each assigned to a different FCC layer so shadows
%   never overlap), the total count scales linearly:
%
%       N_total ≈ n_planes × (π × Ry × Rz) / A_cell
%
%   Inverting for a circular blob (Ry = Rz = R):
%
%       R = s × sqrt( N / (π × n_planes × √3/2) )
%
%   For an elliptical blob with aspect ratio k = Rz/Ry:
%
%       Ry = sqrt( N × A_cell / (π × n_planes × k) )
%       Rz = k × Ry
%
% APPROXIMATION NOTE
%   This formula gives a theoretical maximum (the full ellipse, no clipping).
%   The actual count from place_blob_lattice will typically be LOWER because:
%     1. The Glasgow stability envelope clips some positions at the ellipse edge.
%     2. Hex rows near the ellipse boundary are partially outside.
%   Expect the actual count to be ~10–25% below the analytic estimate,
%   depending on how tightly the blob fits inside the stability envelope.
%   Use this script to choose a starting [Ry, Rz], then run
%   place_blob_lattice (or the full constellation generator) to verify.
%
% USAGE
%   1. Set target_N_values, sail_radius_km, min_buffer_km, and optionally
%      n_planes_range and aspect_ratios in the CONFIGURATION section below.
%   2. Run the script.
%   3. Read the output table and figure to choose your blob dimensions.
%   4. Copy the chosen [Ry, Rz] into input_Constellation_V2_Parameters.m.
%
% AUTHORS
%   Planetary Sunshade Foundation

clear; clc;
fprintf('\n');
fprintf('================================================================\n');
fprintf('  INVERSE LATTICE CALCULATOR: Target count → Ellipse radii\n');
fprintf('================================================================\n\n');


% =========================================================================
% CONFIGURATION — edit these values
% =========================================================================

% Target shade counts to evaluate.
% Can be a single value or a vector for comparison.
target_N_values = [4810];

% Lattice parameters — must match input_Constellation_V2_Parameters.m
sail_radius_km = 20;      % physical radius of each heliogyro (km)
min_buffer_km  = 45;      % edge-to-edge clearance within a plane (km)

% Range of plane counts to tabulate.
n_planes_range = 1:10;

% Ellipse aspect ratios (Rz/Ry) to tabulate alongside the circular case.
% 1.0 = circular, 0.6 = current default (Ry=3000, Rz=1800).
aspect_ratios = [1.0, 0.8, 0.6, 0.5];

% Primary target used for the detailed figure (first entry above).
primary_target = target_N_values(1);


% =========================================================================
% DERIVED GEOMETRY — not user-settable
% =========================================================================

% Hex center-to-center spacing
s = 2 * sail_radius_km + min_buffer_km;

% Area of one hex cell (the region "owned" by each shade)
A_cell = s^2 * sqrt(3) / 2;   % km²

% Maximum number of non-overlapping FCC layers given this spacing.
% Cross-layer minimum distance = s/√3. Must exceed 2r for LOS clearance.
cross_layer_dist = s / sqrt(3);
if cross_layer_dist > 2 * sail_radius_km
    max_layers_feasible = 3;
elseif s > 2 * sail_radius_km
    max_layers_feasible = 2;
else
    max_layers_feasible = 1;
end

fprintf('  Sail radius:          %.0f km\n', sail_radius_km);
fprintf('  Min buffer:           %.0f km (edge-to-edge, within plane)\n', min_buffer_km);
fprintf('  Hex spacing (s):      %.0f km (center-to-center)\n', s);
fprintf('  Hex cell area:        %.1f km²\n', A_cell);
fprintf('  Max FCC layers:       %d  ', max_layers_feasible);
if max_layers_feasible == 3
    fprintf('(cross-layer dist = %.1f km > 2r = %.0f km  ✓ A-B-C)\n', ...
            cross_layer_dist, 2*sail_radius_km);
elseif max_layers_feasible == 2
    fprintf('(s = %.0f km > 2r = %.0f km, but s/√3 = %.1f ≤ 2r — only A-B)\n', ...
            s, 2*sail_radius_km, cross_layer_dist);
else
    fprintf('(s = %.0f km ≤ 2r = %.0f km — single layer only)\n', ...
            s, 2*sail_radius_km);
end
fprintf('\n');
fprintf('  Geometry note: FCC layers multiply the count only if n_planes >\n');
fprintf('  1. The layer capacity (1,2,3) caps how many planes add NEW shades\n');
fprintf('  without line-of-sight overlap. After n_planes > max_layers, the\n');
fprintf('  same Y-Z offsets repeat but shades still stack on distinct X planes.\n');
fprintf('\n');


% =========================================================================
% FUNCTION: analytic shade count estimate
% =========================================================================
%   N_approx(Ry, Rz, np) = np * pi * Ry * Rz / A_cell
%
% Inverted for circular blob:
%   R_circular(N, np) = sqrt(N * A_cell / (pi * np))
%
% Inverted for ellipse with aspect k:
%   Ry(N, np, k) = sqrt(N * A_cell / (pi * np * k))
%   Rz = k * Ry

R_circular = @(N, np)    sqrt(N .* A_cell ./ (pi .* np));
Ry_ellipse = @(N, np, k) sqrt(N .* A_cell ./ (pi .* np .* k));


% =========================================================================
% TABLE: Radii for each target count × n_planes
% =========================================================================

for ti = 1:numel(target_N_values)
    N_target = target_N_values(ti);

    fprintf('----------------------------------------------------------------\n');
    fprintf('  TARGET: %d shades\n', N_target);
    fprintf('----------------------------------------------------------------\n\n');

    % Header
    aspect_labels = arrayfun(@(k) sprintf('k=%.1f', k), aspect_ratios, ...
                             'UniformOutput', false);

    fprintf('  %-9s  %-12s', 'n_planes', 'Circular R');
    for ai = 1:numel(aspect_ratios)
        fprintf('  %-22s', aspect_labels{ai});
    end
    fprintf('\n');

    fprintf('  %-9s  %-12s', '', '(km)');
    for ai = 1:numel(aspect_ratios)
        fprintf('  %-10s  %-10s', 'Ry (km)', 'Rz (km)');
    end
    fprintf('\n');

    fprintf('  %s\n', repmat('-', 1, 9 + 4 + 12 + numel(aspect_ratios) * 24));

    for np = n_planes_range
        R_c = R_circular(N_target, np);
        fprintf('  %-9d  %-12.0f', np, R_c);

        for ai = 1:numel(aspect_ratios)
            k = aspect_ratios(ai);
            Ry = Ry_ellipse(N_target, np, k);
            Rz = k * Ry;
            fprintf('  %-10.0f  %-10.0f', Ry, Rz);
        end
        fprintf('\n');
    end

    fprintf('\n');
end


% =========================================================================
% QUICK VERIFICATION: show analytic vs approx for current default settings
% =========================================================================

fprintf('================================================================\n');
fprintf('  VERIFICATION REFERENCE\n');
fprintf('  Current default in input_Constellation_V2_Parameters.m:\n');
fprintf('    ellipse_radii = [3000, 1800] km, n_planes = 3\n');
fprintf('    → Analytic estimate: %.0f shades\n', ...
        3 * pi * 3000 * 1800 / A_cell);
fprintf('  (Run place_blob_lattice with these settings to get exact count.)\n');
fprintf('================================================================\n\n');


% =========================================================================
% FIGURE 1: Radius vs n_planes for all target counts (circular blob)
% =========================================================================

fig1 = figure('Name', 'Inverse Lattice: Radius vs Planes', ...
              'NumberTitle', 'off', 'Position', [50, 100, 900, 500]);

colors = lines(numel(target_N_values));
hold on;

for ti = 1:numel(target_N_values)
    N_t = target_N_values(ti);
    R_vals = R_circular(N_t, n_planes_range);
    plot(n_planes_range, R_vals, 'o-', ...
         'Color', colors(ti,:), 'LineWidth', 2, 'MarkerSize', 7, ...
         'MarkerFaceColor', colors(ti,:), ...
         'DisplayName', sprintf('%s shades', num2sepstr(N_t)));
end

% Mark the max-layers boundary (after this n_planes, no new LOS geometry)
if max_layers_feasible < max(n_planes_range)
    xline(max_layers_feasible + 0.5, '--k', 'LineWidth', 1.2, ...
          'Label', sprintf('max layers (%d)', max_layers_feasible), ...
          'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off');
end

xlabel('Number of planes (n\_planes)', 'FontSize', 12);
ylabel('Circular blob radius R (km)', 'FontSize', 12);
title(sprintf(['Circular blob radius needed to reach target shade count\n' ...
               '(s = %d km, r = %d km; analytic approximation)'], ...
              s, sail_radius_km), 'FontSize', 11);
legend('Location', 'northeast', 'FontSize', 10);
grid on;
box on;
xlim([min(n_planes_range) - 0.3, max(n_planes_range) + 0.3]);
xticks(n_planes_range);


% =========================================================================
% FIGURE 2: Ry vs Rz for the primary target count, parameterised by
%           n_planes and aspect ratio
% =========================================================================

fig2 = figure('Name', sprintf('Inverse Lattice: Ellipse shapes for N=%d', primary_target), ...
              'NumberTitle', 'off', 'Position', [100, 50, 900, 600]);

hold on;

% Shade a reference region: solar disk at L1 (~10980 km radius)
solar_disk_km = 10980;
fill([0, solar_disk_km, solar_disk_km, 0], ...
     [0, 0, solar_disk_km, solar_disk_km], ...
     [0.95, 0.95, 0.85], 'EdgeColor', 'none', 'DisplayName', 'Solar disk (10980 km)');

n_planes_highlight = [1, 2, 3, 5, 7, 10];
colors2 = lines(numel(n_planes_highlight));

for ni = 1:numel(n_planes_highlight)
    np = n_planes_highlight(ni);
    Ry_vals = [];
    Rz_vals = [];
    for ai = 1:numel(aspect_ratios)
        k  = aspect_ratios(ai);
        Ry = Ry_ellipse(primary_target, np, k);
        Rz = k * Ry;
        Ry_vals(end+1) = Ry; %#ok<AGROW>
        Rz_vals(end+1) = Rz; %#ok<AGROW>
    end
    % Also add the circular point
    R_c = R_circular(primary_target, np);
    Ry_vals(end+1) = R_c;
    Rz_vals(end+1) = R_c;
    [Ry_vals, sort_idx] = sort(Ry_vals);
    Rz_vals = Rz_vals(sort_idx);

    plot(Ry_vals, Rz_vals, 'o-', ...
         'Color', colors2(ni,:), 'LineWidth', 1.8, 'MarkerSize', 7, ...
         'MarkerFaceColor', colors2(ni,:), ...
         'DisplayName', sprintf('%d plane%s', np, 's' * (np > 1)));
end

% Mark the current default
plot(3000, 1800, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'k', ...
     'DisplayName', 'Current default [3000, 1800]');

% Iso-aspect-ratio reference lines
for ai = 1:numel(aspect_ratios)
    k = aspect_ratios(ai);
    R_line = linspace(0, solar_disk_km * 1.1, 100);
    plot(R_line, k * R_line, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, ...
         'HandleVisibility', 'off');
    text(solar_disk_km * 0.85, k * solar_disk_km * 0.85, ...
         sprintf('k=%.1f', k), 'FontSize', 8, 'Color', [0.5 0.5 0.5], ...
         'HorizontalAlignment', 'center');
end

xlabel('Semi-axis R_y (km)', 'FontSize', 12);
ylabel('Semi-axis R_z (km)', 'FontSize', 12);
title(sprintf(['Ellipse shapes giving ≈%s shades (analytic)\n' ...
               'Each line = one n\\_planes; points = different aspect ratios'], ...
              num2sepstr(primary_target)), 'FontSize', 11);
legend('Location', 'northwest', 'FontSize', 10);
grid on; box on; axis equal;
xlim([0, solar_disk_km * 1.1]);
ylim([0, solar_disk_km * 1.1]);


% =========================================================================
% FIGURE 3: Sensitivity — how does count vary around a chosen radius?
% =========================================================================
%   Fixed n_planes=3, circular blob. Sweep R and show N vs R.

fig3 = figure('Name', 'N vs Radius sensitivity', ...
              'NumberTitle', 'off', 'Position', [150, 50, 800, 450]);

R_sweep = linspace(100, 12000, 500);   % km
hold on;

np_sweep_list = [1, 2, 3, 5];
col3 = lines(numel(np_sweep_list));

for ni = 1:numel(np_sweep_list)
    np = np_sweep_list(ni);
    N_analytic = np .* pi .* R_sweep.^2 ./ A_cell;
    plot(R_sweep, N_analytic / 1e3, '-', ...
         'Color', col3(ni,:), 'LineWidth', 2, ...
         'DisplayName', sprintf('%d plane%s', np, 's' * (np > 1)));
end

% Mark each target count with horizontal lines
for ti = 1:numel(target_N_values)
    yline(target_N_values(ti) / 1e3, '--', ...
          'Color', [0.4 0.4 0.4], 'LineWidth', 1.0, ...
          'Label', sprintf('%s', num2sepstr(target_N_values(ti))), ...
          'LabelHorizontalAlignment', 'left', ...
          'HandleVisibility', 'off');
end

xlabel('Circular blob radius R (km)', 'FontSize', 12);
ylabel('Total shades (thousands)', 'FontSize', 12);
title(sprintf(['Analytic shade count vs. blob radius\n' ...
               '(circular blob, s = %d km, r = %d km)'], s, sail_radius_km), ...
      'FontSize', 11);
legend('Location', 'northwest', 'FontSize', 10);
grid on; box on;
xlim([0, max(R_sweep)]);


% =========================================================================
% HELPER: number → comma-separated string  (e.g. 12345 → '12,345')
% =========================================================================

function s_out = num2sepstr(n)
    % Format integer n with thousands separators.
    raw = num2str(n);
    groups = {};
    while numel(raw) > 3
        groups{end+1} = raw(end-2:end); %#ok<AGROW>
        raw = raw(1:end-3);
    end
    groups{end+1} = raw;
    s_out = strjoin(fliplr(groups), ',');
end
