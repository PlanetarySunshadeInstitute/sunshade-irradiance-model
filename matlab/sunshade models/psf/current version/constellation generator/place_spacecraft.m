% place_spacecraft.m
%
% Generates spacecraft positions for a constellation near the Sun-Earth L1
% point, distributed across multiple X-depth planes within the Glasgow
% stability envelope and the user-specified constellation footprint.
%
% USAGE
%   positions = place_spacecraft(params, envelope)
%
% For uniform mode, constellation_radius_km is the primary size constraint
% (circle in Y-Z). For polar mode, craft are placed in per-pole ellipses.
%
% The Glasgow envelope is a secondary constraint — it defines where craft
% CAN go from a propulsion standpoint. At the optimal zone the envelope
% is ~222,000 km in Z, far larger than the useful footprint.
%
% Minimum separation is enforced in 3D. Cross-plane separation is
% guaranteed by plane_spacing_km >> min_buffer_km. Within-plane
% separation is enforced by rejection resampling of violating craft.
%
% Authors: Planetary Sunshade Foundation

function [positions, n_violations] = place_spacecraft(params, envelope)

fprintf('[place] Generating %d spacecraft (%s pattern, %d planes)...\n', ...
        params.N, params.pattern, params.n_planes);

PX          = zeros(params.N, 1);
PY          = zeros(params.N, 1);
PZ          = zeros(params.N, 1);
plane_id    = zeros(params.N, 1);
n_violations = 0;

cursor = 1;

for p = 1:params.n_planes

    X_plane = params.X_planes_km(p);
    N_plane = params.N_per_plane(p);

    % Envelope semi-axes at this X (secondary constraint)
    a_env = interp1(envelope.X_km, envelope.a_km, X_plane, 'linear', 'extrap');
    b_env = interp1(envelope.X_km, envelope.b_km, X_plane, 'linear', 'extrap');
    a_env = max(a_env, 1);
    b_env = max(b_env, 1);

    switch params.pattern
        case 'uniform'
            R = params.constellation_radius_km;
            [Y, Z] = pattern_uniform(N_plane, R, a_env, b_env, ...
                                     envelope.safety_margin, params);
            [Y, Z, n_v] = enforce_min_buffer_uniform(Y, Z, R, ...
                                     params.min_buffer_km, a_env, b_env, ...
                                     envelope.safety_margin, params);
        case 'polar'
            [Y, Z] = pattern_polar(N_plane, a_env, b_env, ...
                                    envelope.safety_margin, params, p);
            [Y, Z, n_v] = enforce_min_buffer_polar(Y, Z, ...
                                     params.min_buffer_km, a_env, b_env, ...
                                     envelope.safety_margin, params);
        otherwise
            error('[place] Unknown pattern: "%s"', params.pattern);
    end

    % Enforce minimum separation within this plane
    n_violations = n_violations + n_v;

    idx           = cursor : cursor + N_plane - 1;
    PX(idx)       = X_plane;
    PY(idx)       = Y;
    PZ(idx)       = Z;
    plane_id(idx) = p;
    cursor        = cursor + N_plane;

end

% Final envelope validation
n_invalid = sum(~check_position_valid(envelope, PY, PZ, PX));
if n_invalid > 0
    warning('[place] %d / %d positions outside Glasgow stability envelope.', ...
            n_invalid, params.N);
else
    fprintf('[place] All positions valid within stability envelope.\n');
end

positions.PX       = PX;
positions.PY       = PY;
positions.PZ       = PZ;
positions.NX       = ones(params.N, 1);
positions.NY       = zeros(params.N, 1);
positions.NZ       = zeros(params.N, 1);
positions.plane_id = plane_id;
positions.pattern  = params.pattern;

fprintf('[place] Done.\n\n');

end


% =========================================================================
%  PATTERN: UNIFORM
% =========================================================================

function [Y, Z] = pattern_uniform(N, R, a_env, b_env, safety_margin, params)
% PATTERN_UNIFORM
%
% Places N craft within the constellation footprint (radius R).
% Profile is either flat ('uniform') or centre-peaked ('gaussian'),
% set by params.footprint_profile.
%
% All positions also checked against the Glasgow ellipse envelope.

Y = zeros(N, 1);
Z = zeros(N, 1);
filled = 0;
batch  = ceil(N * 2);

while filled < N

    if strcmp(params.footprint_profile, 'gaussian')
        % Draw from 2D isotropic Gaussian, sigma = footprint_sigma_km
        sig   = params.footprint_sigma_km;
        Y_c   = sig * randn(batch, 1);
        Z_c   = sig * randn(batch, 1);
    else
        % Uniform in bounding square, reject outside circle
        Y_c = (2 * rand(batch, 1) - 1) * R;
        Z_c = (2 * rand(batch, 1) - 1) * R;
    end

    % Accept if inside footprint circle
    in_footprint = (Y_c.^2 + Z_c.^2) <= R^2;

    % Accept if inside Glasgow envelope ellipse
    in_envelope  = (Y_c / a_env).^2 + (Z_c / b_env).^2 <= safety_margin;

    accept = in_footprint & in_envelope;
    Y_in   = Y_c(accept);
    Z_in   = Z_c(accept);

    n_take = min(N - filled, numel(Y_in));
    Y(filled+1 : filled+n_take) = Y_in(1:n_take);
    Z(filled+1 : filled+n_take) = Z_in(1:n_take);
    filled = filled + n_take;

end

end


% =========================================================================
%  PATTERN: POLAR
% =========================================================================

function [Y, Z] = pattern_polar(N, a_env, b_env, safety_margin, params, plane_idx)
% PATTERN_POLAR
%
% If explicit north/south per-plane counts are provided, use ellipses-only
% and place exact counts for this plane. Otherwise fall back to the legacy
% polar_fraction + hemisphere mix.
use_explicit_polar = isfield(params, 'polar_N_north_per_plane') && ...
                     isfield(params, 'polar_N_south_per_plane') && ...
                     numel(params.polar_N_north_per_plane) >= plane_idx && ...
                     numel(params.polar_N_south_per_plane) >= plane_idx;

if use_explicit_polar
    N_north = params.polar_N_north_per_plane(plane_idx);
    N_south = params.polar_N_south_per_plane(plane_idx);

    [Y_n, Z_n] = sample_polar_ellipse_points(N_north, +params.polar_center_z_km, ...
                                             a_env, b_env, safety_margin, params);
    [Y_s, Z_s] = sample_polar_ellipse_points(N_south, -params.polar_center_z_km, ...
                                             a_env, b_env, safety_margin, params);

    Y = [Y_n; Y_s];
    Z = [Z_n; Z_s];
else
    % Legacy behavior: mix uniform footprint with per-pole ellipse targeting
    % via polar_fraction + hemisphere.
    N_uniform  = round(N * (1 - params.polar_fraction));
    N_targeted = N - N_uniform;

    R = params.constellation_radius_km;
    [Y_u, Z_u] = pattern_uniform(N_uniform, R, a_env, b_env, safety_margin, params);

    switch params.hemisphere
        case 'north'
            N_north = N_targeted;
            N_south = 0;
            [Y_u, Z_u] = keep_sign_and_refill_uniform(Y_u, Z_u, N_uniform, +1, ...
                                    R, a_env, b_env, safety_margin, params);
        case 'south'
            N_north = 0;
            N_south = N_targeted;
            [Y_u, Z_u] = keep_sign_and_refill_uniform(Y_u, Z_u, N_uniform, -1, ...
                                    R, a_env, b_env, safety_margin, params);
        case 'both'
            N_north = ceil(N_targeted / 2);   % odd extra goes north
            N_south = N_targeted - N_north;
    end

    [Y_n, Z_n] = sample_polar_ellipse_points(N_north, +params.polar_center_z_km, ...
                                             a_env, b_env, safety_margin, params);
    [Y_s, Z_s] = sample_polar_ellipse_points(N_south, -params.polar_center_z_km, ...
                                             a_env, b_env, safety_margin, params);

    Y = [Y_u; Y_n; Y_s];
    Z = [Z_u; Z_n; Z_s];
end

% Final size guard — should never trigger, but catches any off-by-one
if numel(Y) ~= N
    error('[pattern_polar] Output size %d does not match requested N=%d.', ...
          numel(Y), N);
end

function [Y_out, Z_out] = keep_sign_and_refill_uniform(Y_in, Z_in, N, sign_keep, ...
                                              R, a_env, b_env, safety_margin, params)
if N == 0
    Y_out = zeros(0, 1);
    Z_out = zeros(0, 1);
    return;
end

if sign_keep > 0
    keep = Z_in >= 0;
else
    keep = Z_in <= 0;
end
Y_out = Y_in(keep);
Z_out = Z_in(keep);
filled = numel(Y_out);

if filled >= N
    Y_out = Y_out(1:N);
    Z_out = Z_out(1:N);
    return;
end

batch = max(64, ceil((N - filled) * 3));
while filled < N
    [Y_c, Z_c] = sample_footprint_batch(batch, R, a_env, b_env, ...
                                        safety_margin, params);
    if sign_keep > 0
        keep = Z_c >= 0;
    else
        keep = Z_c <= 0;
    end
    Y_c = Y_c(keep);
    Z_c = Z_c(keep);
    if isempty(Y_c); continue; end

    n_take = min(N - filled, numel(Y_c));
    Y_out(filled+1:filled+n_take, 1) = Y_c(1:n_take);
    Z_out(filled+1:filled+n_take, 1) = Z_c(1:n_take);
    filled = filled + n_take;
end
end

end

function [Y, Z] = sample_polar_ellipse_points(N, Zc, a_env, b_env, ...
                                              safety_margin, params)
if N == 0
    Y = zeros(0, 1);
    Z = zeros(0, 1);
    return;
end

Y = zeros(N, 1);
Z = zeros(N, 1);
filled = 0;
batch = max(64, ceil(N * 2));

while filled < N
    [Y_c, Z_c] = sample_polar_ellipse_batch(batch, Zc, a_env, b_env, ...
                                            safety_margin, params);
    if isempty(Y_c); continue; end
    n_take = min(N - filled, numel(Y_c));
    Y(filled+1:filled+n_take) = Y_c(1:n_take);
    Z(filled+1:filled+n_take) = Z_c(1:n_take);
    filled = filled + n_take;
end
end

% =========================================================================
%  HELPER: ENFORCE MINIMUM BUFFER (within a single plane, 3D-aware)
% =========================================================================

function [Y, Z, n_remain] = enforce_min_buffer_uniform(Y, Z, R, min_buffer_km, ...
                                               a_env, b_env, safety_margin, ...
                                               params)
% ENFORCE_MIN_BUFFER
%
% Finds craft pairs closer than min_buffer_km and re-draws the offending
% craft within the footprint circle. Iterates until all violations are
% resolved or max_attempts is reached.
%
% Violations persist only due to random re-draw bad luck, not density.
% Increasing max_attempts resolves almost all cases.

max_attempts = 50;    % increased from 10 — space is sparse, just needs luck
N = numel(Y);

for attempt = 1:max_attempts

    dY           = Y - Y';
    dZ           = Z - Z';
    D            = sqrt(dY.^2 + dZ.^2);
    D(1:N+1:end) = Inf;

    violators = find(any(D < min_buffer_km, 2));
    if isempty(violators); n_remain = 0; return; end

    n_v    = numel(violators);
    filled = 0;
    batch  = ceil(n_v * 2);
    Y_new  = zeros(n_v, 1);
    Z_new  = zeros(n_v, 1);

    while filled < n_v
        [Y_c, Z_c] = sample_footprint_batch(batch, R, a_env, b_env, ...
                                            safety_margin, params);
        n_take = min(n_v - filled, numel(Y_c));
        Y_new(filled+1:filled+n_take) = Y_c(1:n_take);
        Z_new(filled+1:filled+n_take) = Z_c(1:n_take);
        filled = filled + n_take;
    end

    Y(violators) = Y_new;
    Z(violators) = Z_new;

end

% Final check — if violations remain after max_attempts, they are
% exported at invalid positions. This is a placement failure, not a
% density problem. The space can hold far more craft at this buffer.
% Increasing max_attempts in enforce_min_buffer will resolve this.
D_final              = sqrt((Y - Y').^2 + (Z - Z').^2);
D_final(1:N+1:end)   = Inf;
still_violating      = any(D_final < min_buffer_km, 2);
n_remain             = sum(still_violating);

if n_remain > 0
    min_sep = min(D_final(D_final < min_buffer_km));
    warning(['[place] PLACEMENT FAILURE: %d craft could not be placed within ' ...
             'the minimum buffer of %.0f km after %d attempts.\n' ...
             '         These craft are too close to a neighbor in the output ' ...
             'file and should not be used.\n' ...
             '         Closest violation: %.1f km (buffer requires %.0f km).\n' ...
             '         Fix: increase max_attempts in enforce_min_buffer, or ' ...
             'increase min_buffer_km / reduce sail_radius_km.'], ...
             n_remain, min_buffer_km, max_attempts, min_sep, min_buffer_km);
end
end

function [Y, Z, n_remain] = enforce_min_buffer_polar(Y, Z, min_buffer_km, ...
                                               a_env, b_env, safety_margin, ...
                                               params)
max_attempts = 50;
N = numel(Y);

for attempt = 1:max_attempts

    dY           = Y - Y';
    dZ           = Z - Z';
    D            = sqrt(dY.^2 + dZ.^2);
    D(1:N+1:end) = Inf;

    violators = find(any(D < min_buffer_km, 2));
    if isempty(violators); n_remain = 0; return; end

    n_v    = numel(violators);
    filled = 0;
    batch  = ceil(n_v * 2);
    Y_new  = zeros(n_v, 1);
    Z_new  = zeros(n_v, 1);

use_explicit = isfield(params, 'polar_N_north_per_plane');

if strcmp(params.hemisphere, 'both') || use_explicit
    % Preserve north/south identity: resample each violator
    % back into whichever ellipse it originally came from.
    idx_n = find(Z(violators) >= 0);
    idx_s = find(Z(violators) <  0);

    [Y_n, Z_n] = sample_polar_ellipse_points(numel(idx_n), +params.polar_center_z_km, ...
                                             a_env, b_env, safety_margin, params);
    [Y_s, Z_s] = sample_polar_ellipse_points(numel(idx_s), -params.polar_center_z_km, ...
                                             a_env, b_env, safety_margin, params);
    Y_new(idx_n) = Y_n;  Z_new(idx_n) = Z_n;
    Y_new(idx_s) = Y_s;  Z_new(idx_s) = Z_s;
else
    % legacy single-hemisphere resampling
    while filled < n_v
        ...
    end
end

    Y(violators) = Y_new;
    Z(violators) = Z_new;
end

D_final              = sqrt((Y - Y').^2 + (Z - Z').^2);
D_final(1:N+1:end)   = Inf;
still_violating      = any(D_final < min_buffer_km, 2);
n_remain             = sum(still_violating);

if n_remain > 0
    min_sep = min(D_final(D_final < min_buffer_km));
    warning(['[place] PLACEMENT FAILURE: %d craft could not be placed within ' ...
             'the minimum buffer of %.0f km after %d attempts.\n' ...
             '         These craft are too close to a neighbor in the output ' ...
             'file and should not be used.\n' ...
             '         Closest violation: %.1f km (buffer requires %.0f km).\n' ...
             '         Fix: increase max_attempts in enforce_min_buffer, or ' ...
             'increase min_buffer_km / reduce sail_radius_km.'], ...
             n_remain, min_buffer_km, max_attempts, min_sep, min_buffer_km);
end
end


% =========================================================================
%  HELPER: SAMPLE FOOTPRINT BATCH
% =========================================================================

function [Y_in, Z_in] = sample_footprint_batch(batch, R, a_env, b_env, ...
                                               safety_margin, params)
if strcmp(params.footprint_profile, 'gaussian')
    sig = params.footprint_sigma_km;
    Y_c = sig * randn(batch, 1);
    Z_c = sig * randn(batch, 1);
else
    Y_c = (2 * rand(batch, 1) - 1) * R;
    Z_c = (2 * rand(batch, 1) - 1) * R;
end

% In polar mode, keep re-sampling on the selected hemisphere so the
% collision-resolution step cannot reintroduce opposite-side points.
if strcmp(params.pattern, 'polar')
    if strcmp(params.hemisphere, 'north')
        keep_h = Z_c >= 0;
        Y_c = Y_c(keep_h);
        Z_c = Z_c(keep_h);
    elseif strcmp(params.hemisphere, 'south')
        keep_h = Z_c <= 0;
        Y_c = Y_c(keep_h);
        Z_c = Z_c(keep_h);
    end
end

in_footprint = (Y_c.^2 + Z_c.^2) <= R^2;
in_envelope  = (Y_c / a_env).^2 + (Z_c / b_env).^2 <= safety_margin;
accept       = in_footprint & in_envelope;

Y_in = Y_c(accept);
Z_in = Z_c(accept);
end

function [Y_in, Z_in] = sample_polar_ellipse_batch(batch, Zc, a_env, b_env, ...
                                                   safety_margin, params)
% Area-preserving disk-to-ellipse sampling.
theta = 2 * pi * rand(batch, 1);
rho   = sqrt(rand(batch, 1));

Y_c = params.polar_radius_y_km * rho .* cos(theta);
Z_c = Zc + params.polar_radius_z_km * rho .* sin(theta);

in_env = (Y_c / a_env).^2 + (Z_c / b_env).^2 <= safety_margin;
Y_in   = Y_c(in_env);
Z_in   = Z_c(in_env);
end

function [Y_in, Z_in] = sample_polar_resample_batch(batch, a_env, b_env, ...
                                                    safety_margin, params)
switch params.hemisphere
    case 'north'
        [Y_in, Z_in] = sample_polar_ellipse_batch(batch, +params.polar_center_z_km, ...
                                                  a_env, b_env, safety_margin, params);
    case 'south'
        [Y_in, Z_in] = sample_polar_ellipse_batch(batch, -params.polar_center_z_km, ...
                                                  a_env, b_env, safety_margin, params);
    case 'both'
        n_north = ceil(batch / 2);
        n_south = batch - n_north;
        [Y_n, Z_n] = sample_polar_ellipse_batch(n_north, +params.polar_center_z_km, ...
                                                a_env, b_env, safety_margin, params);
        [Y_s, Z_s] = sample_polar_ellipse_batch(n_south, -params.polar_center_z_km, ...
                                                a_env, b_env, safety_margin, params);
        Y_in = [Y_n; Y_s];
        Z_in = [Z_n; Z_s];
end
end


% =========================================================================
%  HELPER: CHECK POSITION VALID
% =========================================================================

function valid = check_position_valid(envelope, Y_km, Z_km, X_km)
a_interp    = interp1(envelope.X_km, envelope.a_km, X_km, 'linear', 'extrap');
b_interp    = interp1(envelope.X_km, envelope.b_km, X_km, 'linear', 'extrap');
a_interp    = max(a_interp, 1);
b_interp    = max(b_interp, 1);
ellipse_val = (Y_km ./ a_interp).^2 + (Z_km ./ b_interp).^2;
valid       = ellipse_val <= envelope.safety_margin;
end