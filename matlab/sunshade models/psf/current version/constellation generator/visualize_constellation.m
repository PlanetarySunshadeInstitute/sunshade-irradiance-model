% visualize_constellation.m
%
% Produces a single diagnostic figure for a generated constellation.
% Four panels (2x2):
%
%   [1] Y-Z face-on view       — constellation geometry, coloured by plane
%   [2] Limb darkening heatmap — solar disk coverage weighted by LD effectiveness
%   [3] 3-D scatter            — full spatial distribution coloured by plane
%   [4] Condensed stats
%
% USAGE
%   visualize_constellation(positions, params)
%   visualize_constellation(positions, params, envelope)
%   visualize_constellation(positions, params, envelope, violations)
%
% INPUTS
%   positions   struct with fields PX, PY, PZ (km, relative to L1),
%               and plane_id (integer 1..n_planes)
%   params      struct — see constellation generator for full field list
%   violations  (optional) integer count of buffer violations from
%               place_spacecraft; shown in red in the stats panel if > 0
%
% Authors: Planetary Sunshade Foundation

function visualize_constellation(positions, params, ~, violations)

if nargin < 4
    violations = 0;
end

% -------------------------------------------------------------------------
% SETUP
% -------------------------------------------------------------------------
n_planes     = params.n_planes;
plane_colors = lines(n_planes);
dot_size     = 2;

theta_circle = linspace(0, 2*pi, 300);
if strcmp(params.pattern, 'polar')
    R = max(params.polar_radius_y_km, params.polar_radius_z_km);
else
    R = params.constellation_radius_km;
end
solar_R      = params.solar_disk_radius_km;

% L1-to-Sun distance (km) — used for solar disk projection only
d_L1_sun = 148478496;   % ~1 AU minus L1 distance

fig = figure('Name', 'Constellation Diagnostic', ...
             'NumberTitle', 'off', ...
             'Position', [80 80 1200 900]);

% =========================================================================
% PANEL 1 — Y-Z face-on
% Zoomed to the constellation extent, no solar disk or envelope clutter.
% =========================================================================
ax1 = subplot(2, 2, 1);
hold on;

for p = 1:n_planes
    idx = positions.plane_id == p;
    scatter(positions.PY(idx) / 1e3, positions.PZ(idx) / 1e3, ...
            dot_size, plane_colors(p,:), '.');
end

axis equal; grid on; box on;
xlabel('Y (thousand km)');
ylabel('Z (thousand km)');
title('Face-on (Y–Z)', 'FontWeight', 'bold');

legend_entries = arrayfun(@(p) sprintf('Plane %d  (X = %.0f km)', p, ...
                 params.X_planes_km(p)), 1:n_planes, 'UniformOutput', false);
legend(legend_entries, 'FontSize', 7, 'Location', 'southoutside', 'NumColumns', 2);

% =========================================================================
% PANEL 2 — Solar disk: limb darkening-weighted coverage heatmap
%
% Rendering approach — two layers:
%   Layer 1 (background): the solar disk itself, coloured by the bare limb
%     darkening profile I(theta)/I(0) = 0.61 + 0.39*cos(theta).  This
%     gives the sun a realistic warm-centre/cooler-limb appearance and
%     provides intuitive context for the sail positions.  Space outside
%     the disk is black.
%   Layer 2 (sails): each sail's projected footprint is filled with the
%     local LD value at that position, rendered with a red (high LD,
%     centre) → blue (low LD, limb) colormap.  Sails near the centre are
%     intrinsically more effective and appear red/hot; limb sails appear
%     blue/cool.
%
% Colorbar scale runs from 0.61 (limb, least effective) to 1.00
% (disk centre, most effective) — the true LD intensity values.
% =========================================================================
ax2 = subplot(2, 2, 2);
hold on;

% --- Grid in normalised solar-radii units ---
n_grid = 200;
r_norm = linspace(-1, 1, n_grid);
[Ygrid, Zgrid] = meshgrid(r_norm, r_norm);
Rgrid   = sqrt(Ygrid.^2 + Zgrid.^2);
on_disk = Rgrid <= 1;

% LD value at every grid cell
cos_theta_grid = sqrt(max(1 - Rgrid.^2, 0));
ld_map         = NaN(n_grid, n_grid);
ld_map(on_disk) = 0.61 + 0.39 .* cos_theta_grid(on_disk);

% --- Layer 1: solar disk background ---
% Use a muted yellow-orange gradient (centre bright, limb darker) to
% suggest the photosphere.  Rendered separately so the sail layer sits
% on top with its own colormap.
sun_bg = NaN(n_grid, n_grid, 3);
sun_bg(:,:,1) = min(1.00,  0.90 + 0.10 .* cos_theta_grid) .* on_disk;  % R
sun_bg(:,:,2) = min(1.00,  0.55 + 0.20 .* cos_theta_grid) .* on_disk;  % G
sun_bg(:,:,3) = zeros(n_grid, n_grid) .* on_disk;                        % B
sun_bg(repmat(~on_disk, [1 1 3])) = NaN;

image(r_norm, r_norm, sun_bg, 'AlphaData', on_disk);
set(ax2, 'Color', [0 0 0]);   % space is black

% --- Layer 2: sail coverage, coloured by LD effectiveness ---
% Build a mask of cells covered by at least one sail, and record the
% LD value at each covered cell (invariant of sail count — we want to
% show *where* on the disk shades are sitting and how effective those
% positions are, not how many sails are stacked).
sail_mask = false(n_grid, n_grid);
PX = positions.PX;
PY = positions.PY;
PZ = positions.PZ;
sail_r = params.sail_radius_km;

for k = 1:numel(PX)
    scale  = d_L1_sun / (d_L1_sun + PX(k));
    y_proj = PY(k) * scale / solar_R;
    z_proj = PZ(k) * scale / solar_R;
    r_proj = sail_r * scale / solar_R;

    dist2     = (Ygrid - y_proj).^2 + (Zgrid - z_proj).^2;
    sail_mask = sail_mask | (dist2 <= r_proj^2 & on_disk);
end

% Render covered cells using LD value; uncovered cells transparent
sail_ld = ld_map;
sail_ld(~sail_mask) = NaN;

% --- Custom colormap: red (LD=1.00, centre, most effective) →
%                      yellow (LD=0.61, limb, least effective) ---
% At the limb, sails visually merge with the yellow sun background —
% reinforcing that they contribute little extra.  At the centre, vivid
% red stands out sharply, flagging high-value positions.
n_cm   = 256;
t      = linspace(0, 1, n_cm)';   % 0 = low LD (limb, yellow), 1 = high LD (centre, red)
cmap_r = ones(n_cm, 1);                        % R always 1
cmap_g = 0.85 - 0.85 .* t;                    % G: 0.85 (yellow) → 0 (red)
cmap_b = zeros(n_cm, 1);                       % B always 0
cmap_sail = max(0, min(1, [cmap_r, cmap_g, cmap_b]));

imagesc(r_norm, r_norm, sail_ld, 'AlphaData', sail_mask);
colormap(ax2, cmap_sail);
cb = colorbar;
cb.Label.String = 'I(\theta)/I(0) = a_0 + a_1·cos\theta   [a_0=0.61, a_1=0.39]';
cb.FontSize      = 7;
clim([0.61 1.00]);
cb.Ticks      = [0.61, 0.61+0.39/2, 1.00];
cb.TickLabels = {'a_0 = 0.61  (limb)', ...
                 'a_0 + a_1/2 = 0.805', ...
                 'a_0 + a_1 = 1.00  (centre)'};

% Solar disk outline only (no footprint ring)
plot(cos(theta_circle), sin(theta_circle), 'w-', 'LineWidth', 0.8);

axis equal; axis([-1.1 1.1 -1.1 1.1]); box on;
set(ax2, 'XTick', [-1 -0.5 0 0.5 1], 'YTick', [-1 -0.5 0 0.5 1]);
xlabel('Y  (solar radii)');
ylabel('Z  (solar radii)');
title('Solar disk: limb darkening-weighted coverage', 'FontWeight', 'bold');

% =========================================================================
% PANEL 3 — 3-D scatter
% Full spatial distribution; coloured by plane.
% =========================================================================
ax3 = subplot(2, 2, 3);
hold on;

for p = 1:n_planes
    idx = positions.plane_id == p;
    scatter3(positions.PX(idx) / 1e3, ...
             positions.PY(idx) / 1e3, ...
             positions.PZ(idx) / 1e3, ...
             dot_size, plane_colors(p,:), '.');
end

grid on; box on;
xlabel('X from L1 (thousand km)');
ylabel('Y (thousand km)');
zlabel('Z (thousand km)');
title('3-D distribution', 'FontWeight', 'bold');
view(-35, 25);

legend_entries3 = arrayfun(@(p) sprintf('Plane %d', p), ...
                  1:n_planes, 'UniformOutput', false);
legend(legend_entries3, 'FontSize', 7, 'Location', 'northeast');

% =========================================================================
% PANEL 4 — Condensed stats
% =========================================================================
subplot(2, 2, 4);
axis off;

% --- Derived numbers ---
disk_area = pi * solar_R^2;
foot_area = pi * R^2;
area_pct  = min(foot_area / disk_area * 100, 100);

% Per-sail limb darkening weight (radial position on projected disk)
r_on_disk = sqrt(positions.PY.^2 + positions.PZ.^2) / solar_R;
r_on_disk = min(r_on_disk, 1);
ld_w      = 0.61 + 0.39 .* sqrt(max(1 - r_on_disk.^2, 0));
mean_ld   = mean(ld_w);

% Radial diagnostics in footprint-normalized units (helps compare profiles)
r_norm_fp = sqrt(positions.PY.^2 + positions.PZ.^2) / max(R, eps);
r_norm_fp = min(r_norm_fp, 1);
r_pct     = prctile(r_norm_fp, [50 80 95]);

% Hemisphere diagnostics
n_north = sum(positions.PZ >= 0);
n_south = sum(positions.PZ < 0);
p_north = 100 * n_north / max(params.N, 1);
p_south = 100 * n_south / max(params.N, 1);

% --- Text block ---
L = {};
L{end+1} = '\bfConstellation summary\rm';
L{end+1} = ' ';
if isfield(params, 'output_filename')
    [~, fname, fext] = fileparts(params.output_filename);
    L{end+1} = 'Source';
    L{end+1} = sprintf('  %s', [fname, fext]);
end
L{end+1} = sprintf('Pattern        %s',              params.pattern);
L{end+1} = sprintf('Craft          %d  (%d planes)', params.N, params.n_planes);
L{end+1} = sprintf('Plane spacing  %.0f km',         params.plane_spacing_km);
L{end+1} = ' ';
L{end+1} = '\bfFootprint\rm';
L{end+1} = sprintf('  Profile      %s',              params.footprint_profile);
if strcmp(params.pattern, 'polar')
    L{end+1} = sprintf('  Polar ellipses  Y %.0f km, Z %.0f km', ...
                       params.polar_radius_y_km, params.polar_radius_z_km);
else
    L{end+1} = sprintf('  Radius       %.0f km  (%.0f%% of disk area)', R, area_pct);
end
if strcmp(params.footprint_profile, 'gaussian')
    L{end+1} = sprintf('  Sigma        %.0f km',     params.footprint_sigma_km);
end
if strcmp(params.pattern, 'uniform')
    L{end+1} = sprintf('  Coverage     %.0f%%',          params.footprint_pct_of_disk);
end
L{end+1} = sprintf('  r/R p50,p80,p95  %.2f, %.2f, %.2f', ...
                   r_pct(1), r_pct(2), r_pct(3));
L{end+1} = ' ';
L{end+1} = '\bfLimb darkening effectiveness\rm';
L{end+1} = sprintf('  Mean LD wt   %.3f  (1.000 = disk centre)', mean_ld);
L{end+1} = sprintf('  Sail radius  %.0f km',         params.sail_radius_km);
L{end+1} = sprintf('  North/South  %.1f%% / %.1f%%', p_north, p_south);
if strcmp(params.pattern, 'polar')
    L{end+1} = ' ';
    L{end+1} = '\bfPolar targeting\rm';
    L{end+1} = sprintf('  Fraction     %.0f%%',      params.polar_fraction * 100);
    L{end+1} = sprintf('  Base layer   %d craft (%.1f%%)', ...
                       params.N_uniform, 100 * params.N_uniform / max(params.N, 1));
    L{end+1} = sprintf('  Center Z     %.0f km (%s)', ...
                         params.polar_center_z_km, params.hemisphere);
end
L{end+1} = ' ';
buf_idx   = numel(L) + 1;
L{end+1}  = '';    % placeholder — buffer line rendered separately below
L{end+1} = ' ';
L{end+1} = sprintf('Generated  %s', datestr(now, 'yyyy-mm-dd HH:MM'));

n_lines   = numel(L);
y_spacing = 1 / (n_lines + 1);

for i = 1:n_lines
    if i == buf_idx; continue; end   % rendered separately
    text(0.04, 1 - i * y_spacing, L{i}, ...
         'Units',            'normalized', ...
         'FontSize',         8.5, ...
         'FontName',         'Courier', ...
         'Interpreter',      'tex', ...
         'VerticalAlignment', 'middle');
end

% --- Buffer line — two calls at fixed x positions; violation count red if > 0 ---
buf_y = 1 - buf_idx * y_spacing;
if violations == 0
    viol_str   = '(0 violations)';
    viol_color = [0 0 0];
else
    viol_str   = sprintf('(%d violations)', violations);
    viol_color = [0.85 0 0];
end
text(0.04, buf_y, sprintf('Buffer    %.0f km', params.min_buffer_km), ...
     'Units', 'normalized', 'FontSize', 8.5, 'FontName', 'Courier', ...
     'Interpreter', 'none', 'VerticalAlignment', 'middle');
text(0.55, buf_y, viol_str, ...
     'Units', 'normalized', 'FontSize', 8.5, 'FontName', 'Courier', ...
     'Interpreter', 'none', 'VerticalAlignment', 'middle', ...
     'Color', viol_color);

end