% generate_constellation.m
%
% Master script for the PSF constellation generator.
% Edit the USER SETTINGS section, then run.
%
% Sequence:
%   1. Load pre-computed Glasgow envelope
%   2. Define and validate constellation parameters
%   3. Place spacecraft
%   4. Visualize — user confirms before export
%   5. Export to excel/psf model/
%
% Authors: Planetary Sunshade Foundation

clear; clc;

% =========================================================================
% USER SETTINGS
% =========================================================================

% --- Paths ---------------------------------------------------------------
envelope_file = fullfile('/Users/morgangoodwin/Desktop/PSF/MatLab/matlab/sunshade models/psf/current version/constellation generator', ...
                          'equilibrium_envelope.mat');
excel_folder  = '/Users/morgangoodwin/Desktop/PSF/MatLab/excel/psf model';

% --- Constellation parameters --------------------------------------------
user_params.pattern                 = 'polar';   % 'uniform' or 'polar'
user_params.N                       = 10000;
user_params.n_planes                = 10;
user_params.plane_spacing_km        = 10000;

user_params.constellation_radius_km = 7000;
user_params.footprint_profile       = 'gaussian';
user_params.footprint_sigma_fraction= 0.5;

user_params.sail_radius_km          = 20;
user_params.min_buffer_km           = 50;

% Polar only (ignored for uniform):
user_params.polar_fraction          = 0.99;
user_params.target_latitude_deg     = 70;
user_params.latitude_band_width_deg = 20;
user_params.hemisphere              = 'both';

% =========================================================================
% RUN — no edits needed below this line
% =========================================================================

envelope  = load_envelope(envelope_file);
params    = define_constellation_params(user_params);
[positions, n_violations] = place_spacecraft(params, envelope);

% Pre-compute expected output filename so the stats panel can display it.
params.output_filename = sprintf('constellation_%s_%dcrafts_%s.xlsx', ...
                                  params.pattern, params.N, ...
                                  datestr(now, 'yyyy-mm-dd'));

visualize_constellation(positions, params, envelope, n_violations);

% Pause for user confirmation before writing the Excel file
input_val = input('Constellation looks correct? Export to Excel? (y/n): ', 's');
if ~strcmpi(input_val, 'y')
    fprintf('Export cancelled. Adjust parameters and re-run.\n');
    return;
end

filepath = export_kinematics_xlsx(positions, params, excel_folder);
fprintf('Ready to run irradiance model against:\n  %s\n', filepath);