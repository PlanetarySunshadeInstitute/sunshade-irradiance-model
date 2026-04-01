% visualize_constellation.m
%
% Produces a single diagnostic figure for a generated constellation.
% Four panels: Y-Z face-on view, X-Y side view, X-Z side view, stats.
%
% USAGE
%   visualize_constellation(positions, params, envelope)
%
% Authors: Planetary Sunshade Foundation

function visualize_constellation(positions, params, envelope)

% -------------------------------------------------------------------------
% SETUP
% -------------------------------------------------------------------------
n_planes     = params.n_planes;
plane_colors = lines(n_planes);   % one color per plane, no toolbox needed
dot_size     = 1;

theta_circle = linspace(0, 2*pi, 300);
R            = params.constellation_radius_km;
solar_R      = params.solar_disk_radius_km;

fig = figure('Name', 'Constellation Geometry', ...
             'NumberTitle', 'off', ...
             'Position', [100 100 1200 500]);

% -------------------------------------------------------------------------
% PANEL 1 — Y-Z face-on (as seen from Sun / Earth)
% -------------------------------------------------------------------------
ax1 = subplot(1, 4, 1);
hold on;

for p = 1:n_planes
    idx = positions.plane_id == p;
    scatter(positions.PY(idx) / 1e3, positions.PZ(idx) / 1e3, ...
            dot_size, plane_colors(p,:), '.');
end

% Solar disk and footprint circles
plot(solar_R * cos(theta_circle) / 1e3, ...
     solar_R * sin(theta_circle) / 1e3, ...
     'k--', 'LineWidth', 0.8);
plot(R * cos(theta_circle) / 1e3, ...
     R * sin(theta_circle) / 1e3, ...
     'r-', 'LineWidth', 0.8);

axis equal; grid on; box on;
xlabel('Y (thousand km)');
ylabel('Z (thousand km)');
title('Face-on (Y-Z)');

% Plane legend
legend_entries = arrayfun(@(p) sprintf('Plane %d  X=%.0fkm', p, ...
                 params.X_planes_km(p)), 1:n_planes, ...
                 'UniformOutput', false);
legend_entries{end+1} = 'Solar disk';
legend_entries{end+1} = 'Footprint';
legend(legend_entries, 'FontSize', 7, 'Location', 'southoutside', ...
       'NumColumns', 2);

% -------------------------------------------------------------------------
% PANEL 2 — X-Y side view
% -------------------------------------------------------------------------
ax2 = subplot(1, 4, 2);
hold on;

for p = 1:n_planes
    idx = positions.plane_id == p;
    scatter(positions.PX(idx) / 1e3, positions.PY(idx) / 1e3, ...
            dot_size, plane_colors(p,:), '.');
end

% Mark plane positions
for p = 1:n_planes
    xline(params.X_planes_km(p) / 1e3, '--', ...
          'Color', plane_colors(p,:), 'LineWidth', 0.8);
end

grid on; box on;
xlabel('X from L1 (thousand km)');
ylabel('Y (thousand km)');
title('Side view (X-Y)');

% -------------------------------------------------------------------------
% PANEL 3 — X-Z side view
% -------------------------------------------------------------------------
ax3 = subplot(1, 4, 3);
hold on;

for p = 1:n_planes
    idx = positions.plane_id == p;
    scatter(positions.PX(idx) / 1e3, positions.PZ(idx) / 1e3, ...
            dot_size, plane_colors(p,:), '.');
end

for p = 1:n_planes
    xline(params.X_planes_km(p) / 1e3, '--', ...
          'Color', plane_colors(p,:), 'LineWidth', 0.8);
end

% For polar: mark Z center offset
if strcmp(params.pattern, 'polar')
    switch params.hemisphere
        case {'north', 'both'}
            yline(params.polar_center_z_km / 1e3, 'k:', 'LineWidth', 1);
        case {'south', 'both'}
            yline(-params.polar_center_z_km / 1e3, 'k:', 'LineWidth', 1);
    end
    if strcmp(params.hemisphere, 'both')
        yline( params.polar_center_z_km / 1e3, 'k:', 'LineWidth', 1);
        yline(-params.polar_center_z_km / 1e3, 'k:', 'LineWidth', 1);
    end
end

grid on; box on;
xlabel('X from L1 (thousand km)');
ylabel('Z (thousand km)');
title('Side view (X-Z)');

% -------------------------------------------------------------------------
% PANEL 4 — Stats
% -------------------------------------------------------------------------
subplot(1, 4, 4);
axis off;

% Build stats text
lines_text = {};
lines_text{end+1} = sprintf('\\bfConstellation summary\\rm');
lines_text{end+1} = ' ';
lines_text{end+1} = sprintf('Pattern         %s', params.pattern);
lines_text{end+1} = sprintf('Total craft     %d', params.N);
lines_text{end+1} = sprintf('Planes          %d', params.n_planes);
lines_text{end+1} = sprintf('Plane spacing   %.0f km', params.plane_spacing_km);
lines_text{end+1} = ' ';
lines_text{end+1} = sprintf('Footprint');
lines_text{end+1} = sprintf('  Profile       %s', params.footprint_profile);
lines_text{end+1} = sprintf('  Radius        %.0f km', R);
lines_text{end+1} = sprintf('  Solar disk    %.0f km', solar_R);
lines_text{end+1} = sprintf('  Coverage      %.0f%%', params.footprint_pct_of_disk);
if strcmp(params.footprint_profile, 'gaussian')
    lines_text{end+1} = sprintf('  Sigma         %.0f km', params.footprint_sigma_km);
end
lines_text{end+1} = ' ';
lines_text{end+1} = sprintf('Sail radius     %.0f km', params.sail_radius_km);
lines_text{end+1} = sprintf('Min buffer      %.0f km', params.min_buffer_km);
lines_text{end+1} = ' ';

% Plane breakdown
lines_text{end+1} = sprintf('Craft per plane');
for p = 1:n_planes
    lines_text{end+1} = sprintf('  Plane %d  X=%.0f km  N=%d', ...
                                  p, params.X_planes_km(p), params.N_per_plane(p));
end

% Polar-specific
if strcmp(params.pattern, 'polar')
    lines_text{end+1} = ' ';
    lines_text{end+1} = sprintf('Polar targeting');
    lines_text{end+1} = sprintf('  Fraction      %.0f%%', ...
                                  params.polar_fraction * 100);
    lines_text{end+1} = sprintf('  Center Z      %.1f km (%s)', ...
                                  params.polar_center_z_km, params.hemisphere);
    lines_text{end+1} = sprintf('  Z center      %.1f km', params.polar_center_z_km);
end

% Date
lines_text{end+1} = ' ';
lines_text{end+1} = sprintf('Generated  %s', datestr(now, 'yyyy-mm-dd HH:MM'));

% Render as evenly spaced text
n_lines   = numel(lines_text);
y_spacing = 1 / (n_lines + 1);
for i = 1:n_lines
    text(0.05, 1 - i * y_spacing, lines_text{i}, ...
         'Units', 'normalized', ...
         'FontSize', 8, ...
         'FontName', 'Courier', ...
         'Interpreter', 'tex', ...
         'VerticalAlignment', 'middle');
end

% -------------------------------------------------------------------------
% LINK X AXES of side views so they zoom together
% -------------------------------------------------------------------------
linkaxes([ax2 ax3], 'x');

end