% lattice_shade_count_calculator.m
%
% AREA-TO-RADII CALCULATOR
% ============================================================================
%
% PURPOSE
%   Given a desired total shading area (km²) per cluster, compute the
%   ellipse_radii [Ry, Rz] to enter in input_Constellation_V2_Parameters.m.
%
%   Lattice geometry parameters (sail_radius, min_buffer, n_planes) are read
%   automatically from input_Constellation_V2_Parameters.m — no manual sync
%   required. Only active clusters are processed.
%
% USAGE
%   1. Set target_area_km2 and aspect_ratio below — one entry per active
%      cluster, in the order they appear in input_Constellation_V2_Parameters.
%   2. Run the script.
%   3. Copy the printed ellipse_radii values into input_Constellation_V2_Parameters.m.
%
% NOTE ON ACCURACY
%   The formula gives an analytic upper bound. The actual area placed by
%   the lattice will be ~10–20% lower due to clipping at the ellipse
%   boundary and the Glasgow stability envelope. If your target area is a
%   floor rather than a ceiling, multiply your target_area_km2 values by
%   ~1.15 before running.
%
% AUTHORS
%   Planetary Sunshade Foundation

clear; clc;

% =========================================================================
% USER SETTINGS — edit these values
% =========================================================================

% Target effective blocking area (km²) per active cluster.
% The script accounts for heliogyro geometry (annular blades) and opacity
% internally when computing craft count.
% One entry per active cluster, in definition order.
target_area_km2 = [3258000, 1090000];

% Aspect ratio Rz/Ry per active cluster (1.0 = circular).
% One entry per active cluster, in definition order.
% Default: all circular.
aspect_ratio = [1.0, 1.0];


% =========================================================================
% AUTO-LOAD PARAMS — do not edit below this line
% =========================================================================

params        = input_Constellation_V2_Parameters();
model_params  = input_Model_Parameters___0___Sr();

% --- Lattice geometry ---
r_outer = params.placement.lattice.sail_radius_km;
buf     = params.placement.lattice.min_buffer_km;
np      = params.placement.lattice.n_planes;

% Hex center-to-center spacing and cell area
s      = 2 * r_outer + buf;
A_cell = s^2 * sqrt(3) / 2;   % km²

% --- Heliogyro effective area per craft ---
r_inner = model_params.shade.heliogyro.radii.inner;
r_outer_mp = model_params.shade.heliogyro.radii.outer;
kappa   = model_params.shade.heliogyro.material_irradiance_absorption;

% Consistency check: sail_radius_km must match heliogyro outer radius
if abs(r_outer - r_outer_mp) > 1e-6
    error(['Mismatch: sail_radius_km = %.4g km in input_Constellation_V2_Parameters ' ...
           'but heliogyro outer radius = %.4g km in input_Model_Parameters. ' ...
           'These must be equal. Please reconcile before running.'], ...
           r_outer, r_outer_mp);
end

% Effective blocking area per craft (annular area × opacity)
A_craft = pi * (r_outer^2 - r_inner^2) * kappa;   % km²

% --- Active clusters ---
all_clusters    = params.clusters;
active_mask     = arrayfun(@(c) c.active, all_clusters);
active_clusters = all_clusters(active_mask);
n_active        = numel(active_clusters);


% =========================================================================
% VALIDATE USER INPUT
% =========================================================================

if numel(target_area_km2) ~= n_active
    error(['target_area_km2 has %d entries but there are %d active clusters. ' ...
           'Please provide one value per active cluster.'], ...
           numel(target_area_km2), n_active);
end

if numel(aspect_ratio) ~= n_active
    error(['aspect_ratio has %d entries but there are %d active clusters. ' ...
           'Please provide one value per active cluster.'], ...
           numel(aspect_ratio), n_active);
end

if any(aspect_ratio <= 0) || any(aspect_ratio > 1)
    error('aspect_ratio values must be in the range (0, 1]. (Rz <= Ry)');
end


% =========================================================================
% COMPUTE AND PRINT
% =========================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('  AREA-TO-RADII CALCULATOR\n');
fprintf('============================================================\n');
fprintf('  Lattice params (input_Constellation_V2_Parameters.m):\n');
fprintf('    sail_radius = %g km  |  min_buffer = %g km  |  n_planes = %d\n', ...
        r_outer, buf, np);
fprintf('    hex spacing = %g km  |  cell area  = %.1f km²\n', s, A_cell);
fprintf('  Heliogyro params (input_Model_Parameters.m):\n');
fprintf('    r_inner = %g km  |  r_outer = %g km  |  kappa = %.4f\n', ...
        r_inner, r_outer, kappa);
fprintf('    effective area per craft = %.1f km²\n', A_craft);
fprintf('------------------------------------------------------------\n\n');

for i = 1:n_active
    A_target = target_area_km2(i);
    k        = aspect_ratio(i);
    name     = active_clusters(i).name;

    % Craft count implied by the effective area target
    N = A_target / A_craft;

    % Invert lattice formula for ellipse semi-axes
    Ry = sqrt(N * A_cell / (pi * np * k));
    Rz = k * Ry;

    fprintf('  %s\n', name);
    fprintf('    target area  : %s km²\n',  num2sepstr(round(A_target)));
    fprintf('    aspect ratio : %.2f  (Rz/Ry)\n', k);
    fprintf('    implied N    : %s craft\n', num2sepstr(round(N)));
    fprintf('    --> ellipse_radii = [%.0f, %.0f] km\n', Ry, Rz);
    fprintf('\n');
end

fprintf('------------------------------------------------------------\n');
fprintf('  NOTE: These are analytic upper bounds. Actual shading area\n');
fprintf('  will be ~10–20%% lower due to ellipse-boundary and envelope\n');
fprintf('  clipping. To treat the target as a floor, multiply your\n');
fprintf('  target_area_km2 values by ~1.15 and re-run.\n');
fprintf('============================================================\n\n');


% =========================================================================
% HELPER: integer → comma-separated string  (e.g. 12345 → '12,345')
% =========================================================================

function s_out = num2sepstr(n)
    raw = num2str(n);
    groups = {};
    while numel(raw) > 3
        groups{end+1} = raw(end-2:end); %#ok<AGROW>
        raw = raw(1:end-3);
    end
    groups{end+1} = raw;
    s_out = strjoin(fliplr(groups), ',');
end
