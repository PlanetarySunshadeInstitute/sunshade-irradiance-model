% place_blob.m
%
% Places spacecraft within a single blob's elliptical footprint, distributed
% across multiple X-depth planes, with minimum buffer enforcement.
%
% USAGE
%   [positions_3xN, n_violations] = place_blob(blob, offset_km, envelope)
%
% INPUTS
%   blob        Struct with fields:
%                 center_position  [X, Y, Z] in km (blob reference center)
%                 ellipse_radii    [Ry, Rz] semi-axes in Y-Z plane (km)
%                 n_craft          Total spacecraft count
%                 n_planes         Number of X-depth planes
%                 plane_spacing_km Spacing between planes (km)
%                 min_buffer_km    Minimum center-to-center separation (km)
%
%   offset_km   [1x3] or [3x1] displacement [dX, dY, dZ] in km, added to
%               blob.center_position. Typically the output of the blob's
%               motion_function at a given day. Use [0,0,0] for initial
%               (static) placement.
%
%   envelope    Glasgow stability envelope struct (from load_envelope.m),
%               with fields: X_km, a_km, b_km, safety_margin.
%               Used to constrain Y-Z positions at each X plane.
%
% OUTPUTS
%   positions_3xN   3-by-N matrix of [PX; PY; PZ] positions in km.
%                   Rows: PX (sunward), PY (tangential), PZ (ecliptic north).
%                   N = blob.n_craft.
%
%   n_violations    Number of craft that could not satisfy min_buffer_km
%                   after placement attempts. Should be 0. A warning is
%                   raised (not an error) if > 0, consistent with V1.
%
% NOTES
%   Craft maintain fixed relative positions across time steps — this
%   function is called once per blob (at offset=[0,0,0]) to establish
%   the fixed geometry. The motion loop in generate_constellation_v2.m
%   applies subsequent offsets as a rigid translation.
%
%   Cross-plane separation is guaranteed by plane_spacing_km >> min_buffer_km.
%   Within-plane separation is enforced by rejection resampling (50 attempts).
%
% Authors: Planetary Sunshade Foundation

function [positions_3xN, n_violations] = place_blob(blob, offset_km, envelope)

offset_km = offset_km(:)';   % ensure 1x3

% Effective blob center after applying offset
cx = blob.center_position(1) + offset_km(1);
cy = blob.center_position(2) + offset_km(2);
cz = blob.center_position(3) + offset_km(3);

N        = blob.n_craft;
n_planes = blob.n_planes;

fprintf('[place_blob] Placing %d craft across %d planes (ellipse Ry=%.0f km, Rz=%.0f km)...\n', ...
        N, n_planes, blob.ellipse_radii(1), blob.ellipse_radii(2));

% -------------------------------------------------------------------------
% X-plane positions (symmetric around cx)
% -------------------------------------------------------------------------
half_span  = (n_planes - 1) / 2;
offsets_x  = (-half_span : half_span) * blob.plane_spacing_km;
X_planes   = cx + offsets_x;

% Guard: plane spacing vs. buffer
if blob.plane_spacing_km < blob.min_buffer_km
    warning(['[place_blob] plane_spacing_km (%.0f km) < min_buffer_km (%.0f km). ' ...
             'Cross-plane collisions are possible.'], ...
             blob.plane_spacing_km, blob.min_buffer_km);
end

% -------------------------------------------------------------------------
% N per plane (largest-remainder allocation)
% -------------------------------------------------------------------------
base_n      = floor(N / n_planes);
remainder   = N - base_n * n_planes;
N_per_plane = base_n * ones(1, n_planes);
mid_plane   = ceil(n_planes / 2);
N_per_plane(mid_plane) = N_per_plane(mid_plane) + remainder;

% -------------------------------------------------------------------------
% Allocate output
% -------------------------------------------------------------------------
PX = zeros(N, 1);
PY = zeros(N, 1);
PZ = zeros(N, 1);

n_violations = 0;
cursor = 1;

for p = 1:n_planes

    Xp    = X_planes(p);
    Np    = N_per_plane(p);

    % Envelope semi-axes at this X depth (secondary stability constraint)
    a_env = interp1(envelope.X_km, envelope.a_km, Xp, 'linear', 'extrap');
    b_env = interp1(envelope.X_km, envelope.b_km, Xp, 'linear', 'extrap');
    a_env = max(a_env, 1);
    b_env = max(b_env, 1);

    % Place craft within ellipse, then enforce buffer
    [Y, Z]     = sample_ellipse(Np, cy, cz, blob.ellipse_radii, a_env, b_env, envelope.safety_margin);
    [Y, Z, nv] = enforce_buffer(Y, Z, cy, cz, blob.ellipse_radii, ...
                                 blob.min_buffer_km, a_env, b_env, envelope.safety_margin);

    n_violations = n_violations + nv;

    idx        = cursor : cursor + Np - 1;
    PX(idx)    = Xp;
    PY(idx)    = Y;
    PZ(idx)    = Z;
    cursor     = cursor + Np;

end

% Final envelope validation
n_invalid = sum(~check_in_envelope(envelope, PY, PZ, PX));
if n_invalid > 0
    warning('[place_blob] %d / %d positions fall outside the Glasgow stability envelope.', ...
            n_invalid, N);
end

positions_3xN = [PX, PY, PZ]';   % 3-by-N

fprintf('[place_blob] Done. Violations: %d.\n', n_violations);

end


% =========================================================================
%  HELPER: SAMPLE ELLIPSE
% =========================================================================

function [Y, Z] = sample_ellipse(N, cy, cz, ellipse_radii, a_env, b_env, safety_margin)
% Area-preserving disk-to-ellipse sampling.
% Samples N points within the ellipse centered at (cy, cz) with semi-axes
% ellipse_radii(1) [Y] and ellipse_radii(2) [Z], also constrained by the
% Glasgow envelope ellipse.

Ry = ellipse_radii(1);
Rz = ellipse_radii(2);

Y = zeros(N, 1);
Z = zeros(N, 1);
filled = 0;
batch  = max(64, ceil(N * 2));

while filled < N
    theta = 2 * pi * rand(batch, 1);
    rho   = sqrt(rand(batch, 1));          % uniform area density in disk

    Yc = cy + Ry * rho .* cos(theta);
    Zc = cz + Rz * rho .* sin(theta);

    % Constrain to Glasgow envelope (stability requirement)
    in_env = (Yc / a_env).^2 + (Zc / b_env).^2 <= safety_margin;

    Yc = Yc(in_env);
    Zc = Zc(in_env);

    n_take = min(N - filled, numel(Yc));
    Y(filled+1 : filled+n_take) = Yc(1:n_take);
    Z(filled+1 : filled+n_take) = Zc(1:n_take);
    filled = filled + n_take;
end

end


% =========================================================================
%  HELPER: ENFORCE MINIMUM BUFFER
% =========================================================================

function [Y, Z, n_remain] = enforce_buffer(Y, Z, cy, cz, ellipse_radii, ...
                                            min_buffer_km, a_env, b_env, safety_margin)
% Finds craft pairs closer than min_buffer_km and resamples the violators
% within the blob ellipse. Runs up to max_attempts iterations.
% Reports remaining violations as a warning (consistent with V1 behavior).

max_attempts = 50;
N = numel(Y);

for attempt = 1:max_attempts

    dY             = Y - Y';
    dZ             = Z - Z';
    D              = sqrt(dY.^2 + dZ.^2);
    D(1:N+1:end)   = Inf;                  % zero out self-distances

    violators = find(any(D < min_buffer_km, 2));
    if isempty(violators)
        n_remain = 0;
        return;
    end

    % Resample violators within the blob ellipse
    n_v    = numel(violators);
    Y_new  = zeros(n_v, 1);
    Z_new  = zeros(n_v, 1);
    filled = 0;
    batch  = max(64, ceil(n_v * 2));

    while filled < n_v
        [Yb, Zb] = sample_ellipse(batch, cy, cz, ellipse_radii, a_env, b_env, safety_margin);
        n_take = min(n_v - filled, numel(Yb));
        Y_new(filled+1:filled+n_take) = Yb(1:n_take);
        Z_new(filled+1:filled+n_take) = Zb(1:n_take);
        filled = filled + n_take;
    end

    Y(violators) = Y_new;
    Z(violators) = Z_new;

end

% Final check after all attempts
D_final            = sqrt((Y - Y').^2 + (Z - Z').^2);
D_final(1:N+1:end) = Inf;
still_violating    = any(D_final < min_buffer_km, 2);
n_remain           = sum(still_violating);

if n_remain > 0
    min_sep = min(D_final(D_final < min_buffer_km));
    warning(['[place_blob] PLACEMENT WARNING: %d craft could not satisfy ' ...
             'min_buffer_km = %.0f km after %d attempts.\n' ...
             '             Closest pair: %.1f km apart.\n' ...
             '             Fix: increase ellipse_radii, reduce n_craft, ' ...
             'or increase min_buffer_km.'], ...
             n_remain, min_buffer_km, max_attempts, min_sep);
end

end


% =========================================================================
%  HELPER: CHECK ENVELOPE
% =========================================================================

function valid = check_in_envelope(envelope, Y_km, Z_km, X_km)
a_interp    = interp1(envelope.X_km, envelope.a_km, X_km, 'linear', 'extrap');
b_interp    = interp1(envelope.X_km, envelope.b_km, X_km, 'linear', 'extrap');
a_interp    = max(a_interp, 1);
b_interp    = max(b_interp, 1);
ellipse_val = (Y_km ./ a_interp).^2 + (Z_km ./ b_interp).^2;
valid       = ellipse_val <= envelope.safety_margin;
end
