% define_constellation_params.m
%
% Validates and returns a complete constellation parameter struct.
% Fill in the template below and pass it to generate_constellation.m.
%
% USAGE
%   params = define_constellation_params(user_params)
%
% All fields have defaults — you only need to specify what you want to
% change. Unrecognised fields trigger a warning. Invalid values error.
%
% PATTERN OPTIONS (build incrementally — others added as developed)
%   'uniform'  : even fill across the ellipse envelope
%   'polar'    : two polar ellipses centered at +/- user-specified Z
%
% COORDINATE CONVENTION
%   All positions are L1-centred synodic km (see preprocess_glasgow_data.m).
%   X positive = sunward. Y tangential. Z out-of-plane (ecliptic north).
%
% Authors: Planetary Sunshade Foundation

function params = define_constellation_params(user_params)

% -------------------------------------------------------------------------
% CONSTANTS (geometry — not user-configurable)
% -------------------------------------------------------------------------
% Solar disk projected radius at the craft plane (optimal zone).
% Beyond this radius, sails shade empty sky and contribute zero flux reduction.
% Derived from: R_sun * (D_earth_to_craft / D_sun_to_craft)
% = 696000 * (2360000 / 147238000) = 10,980 km
SOLAR_DISK_RADIUS_KM = 10980;

% -------------------------------------------------------------------------
% DEFAULTS
% -------------------------------------------------------------------------
d = struct();

% --- Core ---
d.pattern              = 'uniform';          
d.N                    = 8000;
d.n_planes             = 5;
d.plane_spacing_km     = 10000;    % spacing between X-depth planes (km)
                                    % planes centred on X_optimal_km (~838,400 km)

% --- Footprint ---
d.constellation_radius_km  = 8784; % radius of constellation in Y-Z plane (km)
                                    % default = 80% of solar disk radius
                                    % max useful = 10,980 km (full solar disk)
                                    % beyond this, sails shade empty sky
d.footprint_profile        = 'gaussian'; % 'uniform'  : even areal density
                                          % 'gaussian' : peaks at center, falls
                                          %   toward limb — more efficient per
                                          %   sail due to limb darkening
d.footprint_sigma_fraction = 0.4;  % gaussian only: sigma as fraction of
                                    % constellation_radius_km. 0.4 puts ~86%
                                    % of craft within the radius.

% --- Collision avoidance ---
p.sail_radius_km           = 20;     % physical radius of each sail (km)
d.min_buffer_km            = 100;     % minimum edge-to-edge clearance (km)

% --- Polar pattern (ignored for 'uniform') ---
d.polar_fraction              = 0.33;
d.polar_center_z_km           = 6200;
d.hemisphere                  = 'both';
d.polar_radius_y_km           = 2000;
d.polar_radius_z_km           = 1000;

% Explicit north/south ellipse counts (new interface; optional)
% If provided (both non-NaN), the polar pattern will place ellipses-only
% craft with exact north/south counts and will derive hemisphere internally.
d.polar_N_north               = NaN;
d.polar_N_south                = NaN;

% -------------------------------------------------------------------------
% KNOWN FIELDS
% -------------------------------------------------------------------------
known_fields = { ...
    'pattern', 'N', 'n_planes', 'plane_spacing_km', ...
    'constellation_radius_km', 'footprint_profile', 'footprint_sigma_fraction', ...
    'sail_radius_km', 'min_buffer_km', ...
    'polar_fraction', 'polar_center_z_km', ...
    'hemisphere', ...
    'polar_radius_y_km', 'polar_radius_z_km' };

known_fields = [known_fields, {'polar_N_north', 'polar_N_south'}];

% -------------------------------------------------------------------------
% MERGE USER INPUT OVER DEFAULTS
% -------------------------------------------------------------------------
params = d;

if nargin == 0 || isempty(user_params)
    fprintf('[params] No input provided — using all defaults.\n');
else
    user_fields = fieldnames(user_params);
    for i = 1:numel(user_fields)
        f = user_fields{i};
        if ismember(f, known_fields)
            params.(f) = user_params.(f);
        else
            warning('[params] Unrecognised field ignored: "%s"', f);
        end
    end
end

% -------------------------------------------------------------------------
% VALIDATE
% -------------------------------------------------------------------------
errors = {};

% pattern
valid_patterns = {'uniform', 'polar'};
if ~ismember(params.pattern, valid_patterns)
    errors{end+1} = sprintf('pattern must be one of: %s', ...
                             strjoin(valid_patterns, ', '));
end

% N
if ~isnumeric(params.N) || params.N < 1 || params.N ~= round(params.N)
    errors{end+1} = 'N must be a positive integer';
end

% n_planes
if ~isnumeric(params.n_planes) || params.n_planes < 1 || ...
        params.n_planes ~= round(params.n_planes)
    errors{end+1} = 'n_planes must be a positive integer';
end

% plane_spacing_km
if ~isnumeric(params.plane_spacing_km) || params.plane_spacing_km <= 0
    errors{end+1} = 'plane_spacing_km must be a positive number';
end

% polar_fraction
if ~isnumeric(params.polar_fraction) || ...
        params.polar_fraction < 0 || params.polar_fraction > 1
    errors{end+1} = 'polar_fraction must be in [0, 1]';
end

% polar_center_z_km
if ~isnumeric(params.polar_center_z_km) || params.polar_center_z_km < 0
    errors{end+1} = 'polar_center_z_km must be a non-negative number';
end

% hemisphere
valid_hemispheres = {'north', 'south', 'both'};
if ~ismember(params.hemisphere, valid_hemispheres)
    errors{end+1} = sprintf('hemisphere must be one of: %s', ...
                             strjoin(valid_hemispheres, ', '));
end

% polar_N_north / polar_N_south (optional explicit counts)
if isfield(params, 'polar_N_north') && ~isnan(params.polar_N_north)
    if ~isnumeric(params.polar_N_north) || params.polar_N_north < 0 || ...
            params.polar_N_north ~= round(params.polar_N_north)
        errors{end+1} = 'polar_N_north must be a non-negative integer (or NaN to ignore)';
    end
end
if isfield(params, 'polar_N_south') && ~isnan(params.polar_N_south)
    if ~isnumeric(params.polar_N_south) || params.polar_N_south < 0 || ...
            params.polar_N_south ~= round(params.polar_N_south)
        errors{end+1} = 'polar_N_south must be a non-negative integer (or NaN to ignore)';
    end
end

% polar ellipse radii
if ~isnumeric(params.polar_radius_y_km) || params.polar_radius_y_km <= 0
    errors{end+1} = 'polar_radius_y_km must be a positive number';
end
if ~isnumeric(params.polar_radius_z_km) || params.polar_radius_z_km <= 0
    errors{end+1} = 'polar_radius_z_km must be a positive number';
end

% Throw all errors at once so the user can fix everything in one pass
if ~isempty(errors)
    msg = sprintf('  - %s\n', errors{:});
    error('[params] Invalid parameters:\n%s', msg);
end

% constellation_radius_km (uniform mode only)
if strcmp(params.pattern, 'uniform')
    if ~isnumeric(params.constellation_radius_km) || ...
            params.constellation_radius_km <= 0
        errors{end+1} = 'constellation_radius_km must be a positive number';
    end
    if params.constellation_radius_km > SOLAR_DISK_RADIUS_KM
        warning(['[params] constellation_radius_km (%.0f km) exceeds solar disk ' ...
                 'radius (%.0f km). Craft beyond the solar disk shade empty sky.'], ...
                 params.constellation_radius_km, SOLAR_DISK_RADIUS_KM);
    end
end

% footprint_profile
valid_profiles = {'uniform', 'gaussian'};
if ~ismember(params.footprint_profile, valid_profiles)
    errors{end+1} = sprintf('footprint_profile must be one of: %s', ...
                             strjoin(valid_profiles, ', '));
end

% footprint_sigma_fraction
if ~isnumeric(params.footprint_sigma_fraction) || ...
        params.footprint_sigma_fraction <= 0 || ...
        params.footprint_sigma_fraction > 1
    errors{end+1} = 'footprint_sigma_fraction must be in (0, 1]';
end

% sail_radius_km
if ~isnumeric(params.sail_radius_km) || params.sail_radius_km <= 0
    errors{end+1} = 'sail_radius_km must be a positive number';
end

if ~isnumeric(params.min_buffer_km) || params.min_buffer_km <= 0
    errors{end+1} = 'min_buffer_km must be a positive number';
end
if params.min_buffer_km < 2 * params.sail_radius_km
    warning(['[params] min_buffer_km (%.0f km) is less than sail diameter ' ...
             '(2 × %.0f km = %.0f km). min_buffer_km is a center-to-center ' ...
             'spacing, so sails will physically overlap by %.0f km.'], ...
             params.min_buffer_km, params.sail_radius_km, ...
             2 * params.sail_radius_km, ...
             2 * params.sail_radius_km - params.min_buffer_km);
end

% -------------------------------------------------------------------------
% DERIVED QUANTITIES
% -------------------------------------------------------------------------
% Compute these once here so downstream modules don't repeat the geometry.

% N per layer (polar only — ignored for uniform)
use_polar_n_counts = isfield(params, 'polar_N_north') && isfield(params, 'polar_N_south') && ...
                     ~isnan(params.polar_N_north) && ~isnan(params.polar_N_south);

if use_polar_n_counts
    params.N_uniform  = 0;
    params.N_targeted = params.polar_N_north + params.polar_N_south;
    params.N          = params.N_targeted;  % override total N to match explicit counts
    if params.polar_N_north > 0 && params.polar_N_south > 0
        params.hemisphere = 'both';
    elseif params.polar_N_north > 0
        params.hemisphere = 'north';
    elseif params.polar_N_south > 0
        params.hemisphere = 'south';
    else
        error('[params] For polar: polar_N_north + polar_N_south must be > 0');
    end
    params.polar_fraction = 1; % ellipses-only when explicit counts are used
else
    params.N_uniform  = round(params.N * (1 - params.polar_fraction));
    params.N_targeted = params.N - params.N_uniform;
end

% Nominal optimal X location for plane centering (km)
X_optimal_km      = 838400;

% X positions of each plane, centred on X_optimal_km
% e.g. 5 planes at 10,000 km spacing -> offsets [-20000 -10000 0 10000 20000]
half_span        = (params.n_planes - 1) / 2;
plane_offsets    = (-half_span : half_span) * params.plane_spacing_km;
params.X_planes_km = X_optimal_km + plane_offsets;

% N per plane (distribute as evenly as possible; remainder goes to middle)
base_per_plane     = floor(params.N / params.n_planes);
remainder          = params.N - base_per_plane * params.n_planes;
params.N_per_plane = base_per_plane * ones(1, params.n_planes);
mid                = ceil(params.n_planes / 2);
params.N_per_plane(mid) = params.N_per_plane(mid) + remainder;

% Polar explicit north/south counts (if provided):
% Allocate requested north/south ellipse craft counts across planes,
% so that the per-plane totals still sum to params.N_per_plane(p).
use_polar_n_counts = isfield(params, 'polar_N_north') && isfield(params, 'polar_N_south') && ...
                     ~isnan(params.polar_N_north) && ~isnan(params.polar_N_south);
if strcmp(params.pattern, 'polar') && use_polar_n_counts
    N_total = params.N;
    if N_total < 1
        error('[params] For polar: total N must be >= 1.');
    else
        Nn_total = params.polar_N_north;
        % Allocate with largest-remainder method to match exact integer sum.
        ideal_north = Nn_total * (params.N_per_plane ./ N_total);
        north_floor = floor(ideal_north);
        remainder_north = Nn_total - sum(north_floor);

        frac = ideal_north - north_floor;
        north_per_plane = north_floor;
        if remainder_north > 0
            [~, sort_idx] = sort(frac, 'descend');
            for ii = 1:remainder_north
                north_per_plane(sort_idx(ii)) = north_per_plane(sort_idx(ii)) + 1;
            end
        elseif remainder_north < 0
            % Shouldn't happen with floor(), but guard anyway.
            [~, sort_idx] = sort(frac, 'ascend');
            for ii = 1:(-remainder_north)
                north_per_plane(sort_idx(ii)) = north_per_plane(sort_idx(ii)) - 1;
            end
        end

        % Sanity checks: non-negative and sums preserved.
        if any(north_per_plane < 0)
            error('[params] For polar: north per-plane allocation produced negative counts.');
        end
        if sum(north_per_plane) ~= Nn_total
            error('[params] For polar: north per-plane allocation does not sum to polar_N_north.');
        end

        params.polar_N_north_per_plane = north_per_plane;
        params.polar_N_south_per_plane = params.N_per_plane - north_per_plane;
    end
end

% Guard: plane spacing must be at least min_buffer for cross-plane safety
if params.plane_spacing_km < params.min_buffer_km
    warning(['[params] plane_spacing_km (%.0f) is less than min_buffer_km (%.0f). ' ...
             'Cross-plane collisions are possible.'], ...
             params.plane_spacing_km, params.min_buffer_km);
end

% Gaussian sigma in km
if strcmp(params.pattern, 'uniform')
    params.footprint_sigma_km = params.constellation_radius_km * ...
                                 params.footprint_sigma_fraction;
else
    params.footprint_sigma_km = NaN;
end

% Solar disk radius — stored for reference and for export header
params.solar_disk_radius_km = SOLAR_DISK_RADIUS_KM;
if strcmp(params.pattern, 'uniform')
    params.footprint_pct_of_disk = 100 * params.constellation_radius_km / ...
                                    SOLAR_DISK_RADIUS_KM;
else
    params.footprint_pct_of_disk = NaN;
end

% -------------------------------------------------------------------------
% SUMMARY
% -------------------------------------------------------------------------
fprintf('\n--- Constellation Parameters ---\n');
fprintf('  Pattern              : %s\n',   params.pattern);
fprintf('  Total spacecraft     : %d\n',   params.N);
fprintf('  Planes               : %d  (spacing: %.0f km)\n', ...
        params.n_planes, params.plane_spacing_km);
fprintf('  X plane positions    : '); fprintf('%.0f  ', params.X_planes_km);
fprintf('km\n');
fprintf('  N per plane          : '); fprintf('%d  ',   params.N_per_plane);
fprintf('\n');
if strcmp(params.pattern, 'uniform')
    fprintf('  Constellation radius : %.0f km  (%.0f%% of solar disk)\n', ...
            params.constellation_radius_km, params.footprint_pct_of_disk);
end
fprintf('  Footprint profile    : %s', params.footprint_profile);
if strcmp(params.footprint_profile, 'gaussian')
    fprintf('  (sigma = %.0f km)', params.footprint_sigma_km);
end
fprintf('\n');
fprintf('  Sail radius          : %.1f km\n', params.sail_radius_km);
fprintf('  Min buffer           : %.1f km  (edge-to-edge)\n', params.min_buffer_km);

if strcmp(params.pattern, 'polar')
    use_polar_n_counts = isfield(params, 'polar_N_north') && isfield(params, 'polar_N_south') && ...
                         ~isnan(params.polar_N_north) && ~isnan(params.polar_N_south);
    if use_polar_n_counts
        fprintf('  Polar north shades   : %d\n', params.polar_N_north);
        fprintf('  Polar south shades   : %d\n', params.polar_N_south);
    end
    fprintf('  Polar fraction       : %.2f  (%d targeted / %d uniform)\n', ...
            params.polar_fraction, params.N_targeted, params.N_uniform);
    fprintf('  Polar center Z       : %.1f km (%s)\n', ...
            params.polar_center_z_km, params.hemisphere);
    fprintf('  Polar ellipse axes   : Y %.0f km, Z %.0f km\n', ...
            params.polar_radius_y_km, params.polar_radius_z_km);
end
fprintf('--------------------------------\n\n');

end