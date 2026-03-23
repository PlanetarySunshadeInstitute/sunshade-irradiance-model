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
%   'polar'    : uniform base layer + Gaussian targeting toward a latitude
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
d.target_latitude_deg         = 60;
d.latitude_band_width_deg     = 20;
d.hemisphere                  = 'both';

% -------------------------------------------------------------------------
% KNOWN FIELDS
% -------------------------------------------------------------------------
known_fields = { ...
    'pattern', 'N', 'n_planes', 'plane_spacing_km', ...
    'constellation_radius_km', 'footprint_profile', 'footprint_sigma_fraction', ...
    'sail_radius_km', 'min_buffer_km', ...
    'polar_fraction', 'target_latitude_deg', ...
    'latitude_band_width_deg', 'hemisphere' };

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

% target_latitude_deg
if ~isnumeric(params.target_latitude_deg) || ...
        abs(params.target_latitude_deg) > 90
    errors{end+1} = 'target_latitude_deg must be in [-90, 90]';
end

% latitude_band_width_deg
if ~isnumeric(params.latitude_band_width_deg) || ...
        params.latitude_band_width_deg <= 0 || ...
        params.latitude_band_width_deg > 90
    errors{end+1} = 'latitude_band_width_deg must be in (0, 90]';
end

% hemisphere
valid_hemispheres = {'north', 'south', 'both'};
if ~ismember(params.hemisphere, valid_hemispheres)
    errors{end+1} = sprintf('hemisphere must be one of: %s', ...
                             strjoin(valid_hemispheres, ', '));
end

% Throw all errors at once so the user can fix everything in one pass
if ~isempty(errors)
    msg = sprintf('  - %s\n', errors{:});
    error('[params] Invalid parameters:\n%s', msg);
end

% constellation_radius_km
if ~isnumeric(params.constellation_radius_km) || ...
        params.constellation_radius_km <= 0
    errors{end+1} = 'constellation_radius_km must be a positive number';
end
if params.constellation_radius_km > SOLAR_DISK_RADIUS_KM
    warning(['[params] constellation_radius_km (%.0f km) exceeds solar disk ' ...
             'radius (%.0f km). Craft beyond the solar disk shade empty sky.'], ...
             params.constellation_radius_km, SOLAR_DISK_RADIUS_KM);
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
params.N_uniform  = round(params.N * (1 - params.polar_fraction));
params.N_targeted = params.N - params.N_uniform;

% Z center offset for target latitude, derived from shadow geometry:
%   Z_craft = R_earth * sin(lat) / (AU / D_craft_from_sun)
% where D_craft_from_sun = AU - (L1_spice + X_optimal)
% Projection factor computed from standard constants — no SPICE needed here
% because we only need the nominal value for placement geometry.
AU_km             = 149597870.7;
L1_spice_km       = 1521600;
X_optimal_km      = 838400;
D_craft_from_sun  = AU_km - (L1_spice_km + X_optimal_km);
projection_factor = AU_km / D_craft_from_sun;   % ~1.016
earth_radius_km   = 6371;

lat_rad           = deg2rad(params.target_latitude_deg);
params.Z_center_km = earth_radius_km * sin(lat_rad) / projection_factor;

% Z sigma: convert latitude band width to km using same projection
band_rad           = deg2rad(params.latitude_band_width_deg);
params.Z_sigma_km  = earth_radius_km * sin(band_rad) / projection_factor;

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

% Guard: plane spacing must be at least min_buffer for cross-plane safety
if params.plane_spacing_km < params.min_buffer_km
    warning(['[params] plane_spacing_km (%.0f) is less than min_buffer_km (%.0f). ' ...
             'Cross-plane collisions are possible.'], ...
             params.plane_spacing_km, params.min_buffer_km);
end

% Gaussian sigma in km
params.footprint_sigma_km = params.constellation_radius_km * ...
                             params.footprint_sigma_fraction;

% Solar disk radius — stored for reference and for export header
params.solar_disk_radius_km = SOLAR_DISK_RADIUS_KM;
params.footprint_pct_of_disk = 100 * params.constellation_radius_km / ...
                                SOLAR_DISK_RADIUS_KM;

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
fprintf('  Constellation radius : %.0f km  (%.0f%% of solar disk)\n', ...
        params.constellation_radius_km, params.footprint_pct_of_disk);
fprintf('  Footprint profile    : %s', params.footprint_profile);
if strcmp(params.footprint_profile, 'gaussian')
    fprintf('  (sigma = %.0f km)', params.footprint_sigma_km);
end
fprintf('\n');
fprintf('  Sail radius          : %.1f km\n', params.sail_radius_km);
fprintf('  Min buffer           : %.1f km  (edge-to-edge)\n', params.min_buffer_km);

if strcmp(params.pattern, 'polar')
    fprintf('  Polar fraction       : %.2f  (%d targeted / %d uniform)\n', ...
            params.polar_fraction, params.N_targeted, params.N_uniform);
    fprintf('  Target latitude      : %.1f deg (%s)\n', ...
            params.target_latitude_deg, params.hemisphere);
    fprintf('  Latitude band width  : %.1f deg (1-sigma)\n', ...
            params.latitude_band_width_deg);
    fprintf('  Z center offset      : %.1f km\n', params.Z_center_km);
    fprintf('  Z sigma              : %.1f km\n', params.Z_sigma_km);
end
fprintf('--------------------------------\n\n');

end