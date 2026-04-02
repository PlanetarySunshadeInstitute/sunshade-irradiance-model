% generate_constellation_v2.m
%
% Time-varying constellation generator (V2).
% Edit input_Constellation_V2_Parameters.m, then run this script.
%
% Sequence:
%   1. Load parameters and Glasgow envelope
%   2. Place craft once per cluster (fixed internal geometry)
%   3. Display 3-panel diagnostic (Z motion, 3-D layout, cluster definitions)
%      — user confirms before export
%   4. Assemble full time-varying position array (3 x N x N_days)
%   5. Save .mat file with naming convention and alongside .png diagnostic
%
% OUTPUT FORMAT
%   A .mat file (v7.3) containing:
%     positions    [3 x N_total x N_days] double  — PX, PY, PZ in km
%     normals      [3 x N_total x N_days] double  — unit normals (sun-facing [1,0,0])
%     params       struct                          — full parameter record
%     timestamps   [1 x N_days] datetime          — one per day (2024-01-01 base)
%     cluster_ids  [1 x N_total] int32            — which cluster each craft belongs to
%
% NAMING CONVENTION
%   constellation_[label]_[N]crafts_[yyyy-mm-dd]_[NNN].mat
%   where [NNN] increments if the file already exists (001, 002, ...).
%   A _preview.png is always saved; a _diagnostic.png is saved post-export.
%
% Authors: Planetary Sunshade Foundation

clear; clc;

% =========================================================================
% LOAD PARAMETERS AND ENVELOPE
% =========================================================================

paths         = config_paths();
envelope_file = fullfile(paths.constellation, 'equilibrium_envelope.mat');
mat_folder    = paths.excel_folder;

params   = input_Constellation_V2_Parameters();
envelope = load_envelope(envelope_file);

N_clusters = numel(params.clusters);
N_days     = 365 + strcmp(params.year_type, 'leap');
method     = params.placement.method;

fprintf('\n=== V2 Constellation Generator ===\n');
fprintf('  Label       : %s\n', params.label);
fprintf('  Year type   : %s (%d days)\n', params.year_type, N_days);
fprintf('  Clusters    : %d defined\n', N_clusters);
fprintf('  Placement   : %s\n', method);
fprintf('==================================\n\n');

% =========================================================================
% PHASE 1: PLACE CRAFT (once per cluster, fixed geometry)
% =========================================================================
% Inactive clusters (active = false) are skipped — their n_craft is set to
% zero and their position array is left empty. N_total is computed after
% this loop once all actual counts are known.

fprintf('--- Phase 1: Placing craft (%s) ---\n', method);

cluster_positions = cell(1, N_clusters);

for c = 1:N_clusters
    cluster = params.clusters(c);

    % --- Skip inactive clusters ---
    if ~cluster.active
        fprintf('Cluster %d (%s): inactive — skipping.\n', c, cluster.name);
        params.clusters(c).n_craft          = 0;
        params.clusters(c).n_planes         = 0;
        params.clusters(c).plane_spacing_km = 0;
        params.clusters(c).min_buffer_km    = 0;
        params.clusters(c).n_layers         = 0;
        cluster_positions{c}                = zeros(3, 0);
        continue;
    end

    switch method

        case 'lattice'
            % Craft count is determined by geometry — user does not set it.
            lp = params.placement.lattice;
            cluster_for_placement = struct( ...
                'center_position',  cluster.center_position, ...
                'ellipse_radii',    cluster.ellipse_radii, ...
                'sail_radius_km',   lp.sail_radius_km, ...
                'min_buffer_km',    lp.min_buffer_km, ...
                'n_planes',         lp.n_planes, ...
                'plane_spacing_km', lp.plane_spacing_km);

            fprintf('Cluster %d (%s): center [%.0f, %.0f, %.0f] km  (lattice, count TBD)\n', ...
                    c, cluster.name, cluster.center_position(1), ...
                    cluster.center_position(2), cluster.center_position(3));

            [pos_3xN, meta] = place_blob_lattice( ...
                cluster_for_placement, envelope, ...
                'los_margin_km', lp.los_margin_km, ...
                'verbose', false, ...
                'verify_los', false);

            % Write actual count and effective placement params back so
            % Phase 2, Phase 3, and the filename all read consistent values.
            params.clusters(c).n_craft          = meta.n_total;
            params.clusters(c).n_planes         = lp.n_planes;
            params.clusters(c).plane_spacing_km = lp.plane_spacing_km;
            params.clusters(c).min_buffer_km    = lp.min_buffer_km;
            params.clusters(c).n_layers         = meta.n_layers;

            fprintf('  → %d shades placed (%d FCC layers, LOS clearance %.1f km)\n', ...
                    meta.n_total, meta.n_layers, meta.min_los_clearance_km);

        case 'random'
            rp = params.placement.random;
            cluster_for_placement = struct( ...
                'center_position',  cluster.center_position, ...
                'ellipse_radii',    cluster.ellipse_radii, ...
                'n_craft',          cluster.n_craft, ...
                'n_planes',         rp.n_planes, ...
                'plane_spacing_km', rp.plane_spacing_km, ...
                'min_buffer_km',    rp.min_buffer_km);

            fprintf('Cluster %d (%s): center [%.0f, %.0f, %.0f] km, %d craft\n', ...
                    c, cluster.name, cluster.center_position(1), ...
                    cluster.center_position(2), cluster.center_position(3), ...
                    cluster.n_craft);

            [pos_3xN, n_viol] = place_blob(cluster_for_placement, [0, 0, 0], envelope);

            if n_viol > 0
                warning('[generate_constellation_v2] Cluster %d (%s): %d craft could not satisfy min_buffer_km.', ...
                        c, cluster.name, n_viol);
            end

            params.clusters(c).n_planes         = rp.n_planes;
            params.clusters(c).plane_spacing_km = rp.plane_spacing_km;
            params.clusters(c).min_buffer_km    = rp.min_buffer_km;
            params.clusters(c).n_layers         = 1;

        otherwise
            error('[generate_constellation_v2] Unknown placement method: ''%s''. Use ''lattice'' or ''random''.', method);
    end

    cluster_positions{c} = pos_3xN;
end

% Compute N_total and build cluster ID array now that all counts are known.
N_total      = sum([params.clusters.n_craft]);
cluster_ids  = zeros(1, N_total, 'int32');
cursor       = 1;
for c = 1:N_clusters
    n_c = params.clusters(c).n_craft;
    cluster_ids(cursor : cursor+n_c-1) = int32(c);
    cursor = cursor + n_c;
end

fprintf('\nTotal craft placed: %d\n', N_total);

% =========================================================================
% PHASE 2: DIAGNOSTIC FIGURE (3 panels)
% =========================================================================

fprintf('\n--- Phase 2: Diagnostic figure ---\n');

month_starts = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335];
month_names  = {'Jan','Feb','Mar','Apr','May','Jun', ...
                'Jul','Aug','Sep','Oct','Nov','Dec'};

cluster_colors  = lines(N_clusters);
active_clusters = find([params.clusters.n_craft] > 0);

fig_diag = figure('Name', 'V2 Constellation Diagnostic', 'NumberTitle', 'off', ...
                  'Color', [0.12 0.12 0.12], 'Position', [50, 50, 1600, 720]);

tl = tiledlayout(fig_diag, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Diagnostic — %s  |  %d craft  |  %d days', ...
      params.label, N_total, N_days), ...
      'Color', [0.95 0.95 0.95], 'FontSize', 13, 'FontWeight', 'bold');

% -------------------------------------------------------------------------
% Panel 1 (left column, spans 2 rows): Z-axis motion over the year
% -------------------------------------------------------------------------

ax1 = nexttile(tl, 1, [2, 1]);
hold(ax1, 'on');
set(ax1, 'Color', [0.15 0.15 0.15], ...
         'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], ...
         'GridColor', [0.3 0.3 0.3]);
grid(ax1, 'on');

t_fine = linspace(1, N_days, N_days * 4)';

for c = active_clusters
    cluster = params.clusters(c);
    Z_traj  = zeros(numel(t_fine), 1);
    for ki = 1:numel(t_fine)
        off        = cluster.motion_function(t_fine(ki));
        off        = off(:);
        Z_traj(ki) = cluster.center_position(3) + off(3);
    end
    plot(ax1, t_fine, Z_traj, '-', 'Color', cluster_colors(c,:), 'LineWidth', 2, ...
         'DisplayName', sprintf('%s  (N=%d)', cluster.name, cluster.n_craft));
end

ylabel(ax1, 'Z position (km)', 'Color', [0.8 0.8 0.8]);
title(ax1, 'Z-axis motion over year', 'Color', [0.9 0.9 0.9], 'FontWeight', 'bold');
xticks(ax1, month_starts);
xticklabels(ax1, month_names);
xlim(ax1, [1, N_days]);
if numel(active_clusters) > 1
    legend(ax1, 'TextColor', [0.9 0.9 0.9], 'Color', [0.2 0.2 0.2], 'Location', 'best');
end

% -------------------------------------------------------------------------
% Panel 2 (top-right): 3-D spatial layout — Jan 1
% -------------------------------------------------------------------------

ax2 = nexttile(tl, 2);
hold(ax2, 'on');
set(ax2, 'Color', [0.15 0.15 0.15], ...
         'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], 'ZColor', [0.8 0.8 0.8], ...
         'GridColor', [0.3 0.3 0.3]);
grid(ax2, 'on');

ref_day = 1;
for c = active_clusters
    cluster = params.clusters(c);
    off     = cluster.motion_function(ref_day);
    off     = off(:)';
    pos     = cluster_positions{c};
    scatter3(ax2, pos(1,:) + off(1), pos(2,:) + off(2), pos(3,:) + off(3), ...
             1, cluster_colors(c,:), '.', 'DisplayName', cluster.name);
end

xlabel(ax2, 'X (km)', 'Color', [0.8 0.8 0.8]);
ylabel(ax2, 'Y (km)', 'Color', [0.8 0.8 0.8]);
zlabel(ax2, 'Z (km)', 'Color', [0.8 0.8 0.8]);
title(ax2, 'Clusters in 3-D — Jan 1', 'Color', [0.9 0.9 0.9], 'FontWeight', 'bold');
view(ax2, 35, 20);
if numel(active_clusters) > 1
    legend(ax2, 'TextColor', [0.9 0.9 0.9], 'Color', [0.2 0.2 0.2], 'Location', 'best');
end

% -------------------------------------------------------------------------
% Panel 3 (bottom-right): Cluster definitions printout
% -------------------------------------------------------------------------

ax3 = nexttile(tl, 4);
axis(ax3, 'off');
set(ax3, 'Color', [0.15 0.15 0.15]);

txt_lines = {};
txt_lines{end+1} = 'CLUSTER DEFINITIONS';
txt_lines{end+1} = repmat('-', 1, 44);

for c = 1:N_clusters
    cluster = params.clusters(c);
    if ~cluster.active
        txt_lines{end+1} = sprintf('%s  [inactive]', cluster.name);
        txt_lines{end+1} = '';
        continue;
    end
    txt_lines{end+1} = cluster.name;
    txt_lines{end+1} = sprintf('  center    [%.0f, %.0f, %.0f] km', ...
                                cluster.center_position(1), cluster.center_position(2), ...
                                cluster.center_position(3));
    txt_lines{end+1} = sprintf('  ellipse   Ry=%.0f km, Rz=%.0f km', ...
                                cluster.ellipse_radii(1), cluster.ellipse_radii(2));
    if strcmp(method, 'lattice')
        txt_lines{end+1} = sprintf('  craft     %d   across %d planes (%d FCC layers)', ...
                                    cluster.n_craft, cluster.n_planes, cluster.n_layers);
    else
        txt_lines{end+1} = sprintf('  craft     %d   across %d planes', ...
                                    cluster.n_craft, cluster.n_planes);
    end
    txt_lines{end+1} = sprintf('  spacing   %.0f km,  buffer %.0f km', ...
                                cluster.plane_spacing_km, cluster.min_buffer_km);
    txt_lines{end+1} = sprintf('  method    %s', method);
    txt_lines{end+1} = sprintf('  motion    %s', func2str(cluster.motion_function));
    txt_lines{end+1} = '';
end

text(ax3, 0.04, 0.96, txt_lines, ...
     'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'HorizontalAlignment', 'left', ...
     'Color', [0.85 0.85 0.85], 'FontName', 'Courier', 'FontSize', 7.5, ...
     'Interpreter', 'none');
title(ax3, 'Cluster definitions', 'Color', [0.9 0.9 0.9], 'FontWeight', 'bold');

% Save preview now, while the figure handle is guaranteed valid.
% The file lands in mat_folder so you can open it from Finder during review.
preview_base = sprintf('constellation_%s_%dcrafts_%s_preview.png', ...
                       params.label, N_total, datestr(now, 'yyyy-mm-dd'));
preview_path = fullfile(mat_folder, preview_base);
if ~exist(mat_folder, 'dir'),  mkdir(mat_folder);  end
drawnow;
try
    print(fig_diag, preview_path, '-dpng', '-r200');
    fprintf('\nDiagnostic saved: %s\n', preview_path);
catch err
    fprintf('\n(Could not save preview: %s)\n', err.message);
end

% Pause for confirmation before proceeding with export
fprintf('Review the diagnostic figure (and the PNG above) before continuing.\n');
input_val = input('Proceed with full .mat export? (y/n): ', 's');
if ~strcmpi(input_val, 'y')
    fprintf('Export cancelled. Adjust parameters in input_Constellation_V2_Parameters.m and re-run.\n');
    return;
end

% =========================================================================
% PHASE 3: ASSEMBLE TIME-VARYING ARRAY
% =========================================================================

fprintf('\n--- Phase 3: Assembling %d-day position array (%d craft) ---\n', ...
        N_days, N_total);

positions  = zeros(3, N_total, N_days);
normals    = zeros(3, N_total, N_days);
sun_normal = [1; 0; 0];

for d = 1:N_days
    craft_idx = 1;
    for c = 1:N_clusters
        cluster = params.clusters(c);
        offset  = cluster.motion_function(d);
        offset  = offset(:);

        n_c   = cluster.n_craft;
        pos_d = cluster_positions{c};
        pos_d(1,:) = pos_d(1,:) + offset(1);
        pos_d(2,:) = pos_d(2,:) + offset(2);
        pos_d(3,:) = pos_d(3,:) + offset(3);

        positions(:, craft_idx:craft_idx+n_c-1, d) = pos_d;
        normals(:,   craft_idx:craft_idx+n_c-1, d) = repmat(sun_normal, 1, n_c);
        craft_idx = craft_idx + n_c;
    end

    if mod(d, 50) == 0 || d == N_days
        fprintf('  Day %d / %d\n', d, N_days);
    end
end

timestamps = datetime(2024, 1, 1) + days(0 : N_days-1);

% =========================================================================
% PHASE 4: SAVE .MAT FILE
% =========================================================================

fprintf('\n--- Phase 4: Saving .mat file ---\n');

if ~exist(mat_folder, 'dir')
    mkdir(mat_folder);
    fprintf('Created output folder: %s\n', mat_folder);
end

base_name   = sprintf('constellation_%s_%dcrafts_%s', ...
                       params.label, N_total, datestr(now, 'yyyy-mm-dd'));
output_name = resolve_unique_mat_filename(mat_folder, [base_name, '.mat']);
output_path = fullfile(mat_folder, output_name);

save_data.positions    = positions;
save_data.normals      = normals;
save_data.params       = params;
save_data.timestamps   = timestamps;
save_data.cluster_ids  = cluster_ids;

save(output_path, '-struct', 'save_data', '-v7.3');
fprintf('Saved: %s\n', output_path);

fprintf('\nDone. Output: %s\n', output_name);
fprintf('To run the irradiance model, update location_Heliogyro_Kinematics_Data_V2___E___Sr.m\n');
fprintf('to point at: %s\n', output_name);

% =========================================================================
%  HELPER: RESOLVE UNIQUE .MAT FILENAME
% =========================================================================

function unique_name = resolve_unique_mat_filename(folder, desired_name)
[~, base_name, ~] = fileparts(desired_name);
candidate = [base_name, '.mat'];
serial    = 1;
while isfile(fullfile(folder, candidate))
    candidate = sprintf('%s_%03d.mat', base_name, serial);
    serial    = serial + 1;
end
unique_name = candidate;
end
