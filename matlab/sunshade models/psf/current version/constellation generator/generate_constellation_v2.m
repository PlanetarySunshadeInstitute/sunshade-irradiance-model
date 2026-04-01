% generate_constellation_v2.m
%
% Time-varying constellation generator (V2).
% Edit input_Constellation_V2_Parameters.m, then run this script.
%
% Sequence:
%   1. Load parameters and Glasgow envelope
%   2. Place craft once per blob (fixed internal geometry)
%   3. Display 3-panel diagnostic (Z motion, 3-D layout, blob definitions) — user confirms before export
%   4. Assemble full time-varying position array (3 x N x N_days)
%   5. Save .mat file with naming convention and alongside .png diagnostic
%
% OUTPUT FORMAT
%   A .mat file (v7.3) containing:
%     positions   [3 x N_total x N_days] double  — PX, PY, PZ in km
%     normals     [3 x N_total x N_days] double  — unit normals (sun-facing [1,0,0])
%     params      struct                          — full parameter record
%     timestamps  [1 x N_days] datetime          — one per day (2024-01-01 base)
%     blob_ids    [1 x N_total] int32            — which blob each craft belongs to
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

envelope_file = fullfile('/Users/morgangoodwin/Desktop/PSF/MatLab/matlab/sunshade models/psf/current version/constellation generator', ...
                          'equilibrium_envelope.mat');
mat_folder    = '/Users/morgangoodwin/Desktop/PSF/MatLab/excel/psf model';

params   = input_Constellation_V2_Parameters();
envelope = load_envelope(envelope_file);

N_blobs  = numel(params.blobs);
N_days   = 365 + strcmp(params.year_type, 'leap');
N_total  = sum([params.blobs.n_craft]);

fprintf('\n=== V2 Constellation Generator ===\n');
fprintf('  Label       : %s\n', params.label);
fprintf('  Year type   : %s (%d days)\n', params.year_type, N_days);
fprintf('  Blobs       : %d\n', N_blobs);
fprintf('  Total craft : %d\n', N_total);
fprintf('==================================\n\n');

% =========================================================================
% PHASE 1: PLACE CRAFT (once per blob, fixed geometry)
% =========================================================================

fprintf('--- Phase 1: Placing craft ---\n');

blob_positions = cell(1, N_blobs);   % each cell: 3 x blob.n_craft
blob_ids_all   = zeros(1, N_total, 'int32');
cursor         = 1;

for b = 1:N_blobs
    blob = params.blobs(b);
    fprintf('Blob %d/%d: center [%.0f, %.0f, %.0f] km, %d craft\n', ...
            b, N_blobs, blob.center_position(1), blob.center_position(2), ...
            blob.center_position(3), blob.n_craft);

    [pos_3xN, n_viol] = place_blob(blob, [0, 0, 0], envelope);
    blob_positions{b}  = pos_3xN;    % 3 x n_craft (positions at blob center)

    n_b = blob.n_craft;
    blob_ids_all(cursor : cursor+n_b-1) = int32(b);
    cursor = cursor + n_b;
end

% =========================================================================
% PHASE 2: DIAGNOSTIC FIGURE (3 panels)
% =========================================================================

fprintf('\n--- Phase 2: Diagnostic figure ---\n');

% Month tick positions (day-of-year, 1-indexed, normal year)
month_starts = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335];
month_names  = {'Jan','Feb','Mar','Apr','May','Jun', ...
                'Jul','Aug','Sep','Oct','Nov','Dec'};

blob_colors  = lines(N_blobs);
active_blobs = find([params.blobs.n_craft] > 0);

fig_diag = figure('Name', 'V2 Constellation Diagnostic', 'NumberTitle', 'off', ...
                  'Color', [0.12 0.12 0.12], 'Position', [50, 50, 1600, 720]);

tl = tiledlayout(fig_diag, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Diagnostic — %s  |  %d craft  |  %d days', ...
      params.label, N_total, N_days), ...
      'Color', [0.95 0.95 0.95], 'FontSize', 13, 'FontWeight', 'bold');

% -------------------------------------------------------------------------
% Panel 1 (left column, spans 2 rows): Z-axis motion over the year
%   X axis = time (month labels);  Y axis = Z position (km)
% -------------------------------------------------------------------------

ax1 = nexttile(tl, 1, [2, 1]);
hold(ax1, 'on');
set(ax1, 'Color', [0.15 0.15 0.15], ...
         'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], ...
         'GridColor', [0.3 0.3 0.3]);
grid(ax1, 'on');

% Evaluate motion function at fine resolution for a smooth sin curve
t_fine = linspace(1, N_days, N_days * 4)';   % ~1460 points per blob

for b = active_blobs
    blob   = params.blobs(b);
    Z_traj = zeros(numel(t_fine), 1);
    for ki = 1:numel(t_fine)
        off        = blob.motion_function(t_fine(ki));
        off        = off(:);          % ensure column vector
        Z_traj(ki) = blob.center_position(3) + off(3);
    end
    plot(ax1, t_fine, Z_traj, '-', 'Color', blob_colors(b,:), 'LineWidth', 2, ...
         'DisplayName', sprintf('Blob %d  (N=%d)', b, blob.n_craft));
end

ylabel(ax1, 'Z position (km)', 'Color', [0.8 0.8 0.8]);
title(ax1, 'Z-axis motion over year', 'Color', [0.9 0.9 0.9], 'FontWeight', 'bold');
xticks(ax1, month_starts);
xticklabels(ax1, month_names);
xlim(ax1, [1, N_days]);
if numel(active_blobs) > 1
    legend(ax1, 'TextColor', [0.9 0.9 0.9], 'Color', [0.2 0.2 0.2], 'Location', 'best');
end

% -------------------------------------------------------------------------
% Panel 2 (top-right): 3-D spatial layout — planes visible along X axis
% -------------------------------------------------------------------------

ax2 = nexttile(tl, 2);
hold(ax2, 'on');
set(ax2, 'Color', [0.15 0.15 0.15], ...
         'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], 'ZColor', [0.8 0.8 0.8], ...
         'GridColor', [0.3 0.3 0.3]);
grid(ax2, 'on');

ref_day = 1;
for b = active_blobs
    blob = params.blobs(b);
    off  = blob.motion_function(ref_day);
    off  = off(:)';
    pos  = blob_positions{b};   % 3-by-n_craft  [PX; PY; PZ]
    scatter3(ax2, pos(1,:) + off(1), pos(2,:) + off(2), pos(3,:) + off(3), ...
             1, blob_colors(b,:), '.', 'DisplayName', sprintf('Blob %d', b));
end

xlabel(ax2, 'X (km)', 'Color', [0.8 0.8 0.8]);
ylabel(ax2, 'Y (km)', 'Color', [0.8 0.8 0.8]);
zlabel(ax2, 'Z (km)', 'Color', [0.8 0.8 0.8]);
title(ax2, 'Planes in 3-D — Jan 1', 'Color', [0.9 0.9 0.9], 'FontWeight', 'bold');
view(ax2, 35, 20);
if numel(active_blobs) > 1
    legend(ax2, 'TextColor', [0.9 0.9 0.9], 'Color', [0.2 0.2 0.2], 'Location', 'best');
end

% -------------------------------------------------------------------------
% Panel 3 (bottom-right): Blob definitions printout (n_craft = 0 omitted)
% -------------------------------------------------------------------------

ax3 = nexttile(tl, 4);
axis(ax3, 'off');
set(ax3, 'Color', [0.15 0.15 0.15]);

txt_lines = {};
txt_lines{end+1} = 'BLOB DEFINITIONS';
txt_lines{end+1} = repmat('-', 1, 44);

for b = 1:N_blobs
    blob = params.blobs(b);
    if blob.n_craft == 0
        continue;   % omit zero-craft blobs to reduce confusion
    end
    txt_lines{end+1} = sprintf('Blob %d', b);
    txt_lines{end+1} = sprintf('  center    [%.0f, %.0f, %.0f] km', ...
                                blob.center_position(1), blob.center_position(2), ...
                                blob.center_position(3));
    txt_lines{end+1} = sprintf('  ellipse   Ry=%.0f km, Rz=%.0f km', ...
                                blob.ellipse_radii(1), blob.ellipse_radii(2));
    txt_lines{end+1} = sprintf('  craft     %d   across %d planes', ...
                                blob.n_craft, blob.n_planes);
    txt_lines{end+1} = sprintf('  spacing   %.0f km,  buffer %.0f km', ...
                                blob.plane_spacing_km, blob.min_buffer_km);
    txt_lines{end+1} = sprintf('  motion    %s', func2str(blob.motion_function));
    txt_lines{end+1} = '';
end

text(ax3, 0.04, 0.96, txt_lines, ...
     'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'HorizontalAlignment', 'left', ...
     'Color', [0.85 0.85 0.85], 'FontName', 'Courier', 'FontSize', 7.5, ...
     'Interpreter', 'none');
title(ax3, 'Blob definitions  (n_{craft} = 0 omitted)', ...
      'Color', [0.9 0.9 0.9], 'FontWeight', 'bold');

% Pause for confirmation before proceeding with export
fprintf('\nReview the diagnostic figure.\n');
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

positions = zeros(3, N_total, N_days);
normals   = zeros(3, N_total, N_days);

sun_normal = [1; 0; 0];   % all craft sun-facing

for d = 1:N_days
    craft_idx = 1;
    for b = 1:N_blobs
        blob   = params.blobs(b);
        offset = blob.motion_function(d);
        offset = offset(:);           % ensure 3x1

        n_b    = blob.n_craft;
        % Translate fixed geometry by motion offset at day d
        pos_d  = blob_positions{b};
        pos_d(1,:) = pos_d(1,:) + offset(1);
        pos_d(2,:) = pos_d(2,:) + offset(2);
        pos_d(3,:) = pos_d(3,:) + offset(3);

        positions(:, craft_idx:craft_idx+n_b-1, d) = pos_d;
        normals(:,   craft_idx:craft_idx+n_b-1, d) = repmat(sun_normal, 1, n_b);
        craft_idx = craft_idx + n_b;
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

% Ensure output folder exists (creates it if needed, no-op if it does)
if ~exist(mat_folder, 'dir')
    mkdir(mat_folder);
    fprintf('Created output folder: %s\n', mat_folder);
end

% Build filename following V1 convention:
%   constellation_[label]_[N]crafts_[date]_[NNN].mat
base_name    = sprintf('constellation_%s_%dcrafts_%s', ...
                       params.label, N_total, datestr(now, 'yyyy-mm-dd'));
output_name  = resolve_unique_mat_filename(mat_folder, [base_name, '.mat']);
output_path  = fullfile(mat_folder, output_name);

% Struct to save
save_data.positions  = positions;
save_data.normals    = normals;
save_data.params     = params;
save_data.timestamps = timestamps;
save_data.blob_ids   = blob_ids_all;

save(output_path, '-struct', 'save_data', '-v7.3');
fprintf('Saved: %s\n', output_path);

% Save preview figure as _preview.png alongside the .mat
[folder, base, ~] = fileparts(output_path);
preview_path = fullfile(folder, [base, '_preview.png']);
drawnow;   % flush renderer before capture
try
    print(fig_diag, preview_path, '-dpng', '-r200');
    fprintf('Saved preview: %s\n', preview_path);
catch err
    warning('Could not save preview image: %s\n  %s', preview_path, err.message);
end

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
