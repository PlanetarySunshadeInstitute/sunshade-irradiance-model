% preprocess_glasgow_data.m
%
% ONE-TIME preprocessing script. Reads the raw Glasgow equilibrium region
% Excel file, converts from Sun-origin AU synodic coordinates to L1-centred
% km, fits the ellipse envelope, and saves equilibrium_envelope.mat.
%
% Run this script:
%   - Once when you first set up the repository
%   - Again if Glasgow provides an updated data file
%
% Output: equilibrium_envelope.mat
%   Saved to the same folder as this script. Commit it to the repo so
%   other users never need the raw Excel file to run the generator.
%
% -------------------------------------------------------------------------
% COORDINATE SYSTEM CONTRACT
% -------------------------------------------------------------------------
% This script and the constellation generator use an L1-centred synodic
% frame (Option 3 in the coordinate system analysis):
%
%   Origin : Sun-Earth L1 point
%   X      : positive sunward (toward Sun, away from Earth)
%   Y      : tangential (orbital plane, direction of Earth's motion)
%   Z      : out-of-plane (ecliptic north)
%   Units  : kilometres
%
% Craft positions in this frame are stored in the Excel kinematics file
% as small offsets from L1 (typically <1M km). They are treated as fixed
% — the constellation geometry does not change with the season.
%
% The irradiance model is responsible for adding L1(t) from SPICE at
% runtime to recover absolute inertial positions. This is where
% aphelion/perihelion and annual solar radiation pressure variation
% properly enter the calculation. The constellation generator does not
% need to know about Earth's orbital eccentricity.
%
% TWO L1 VALUES ARE STORED IN THE OUTPUT:
%   L1_glasgow_km : Hill sphere approximation — almost certainly what
%                   Glasgow used to generate their equilibrium data.
%                   Used here to convert Glasgow AU coordinates to km
%                   offsets, ensuring the envelope geometry is
%                   self-consistent with the source data.
%
%   L1_spice_km   : SPICE-validated value from the PSF irradiance model
%                   (distance_PL1___Sr___Sr.m). Stored for reference
%                   and used only to compute X_optimal_km, since that
%                   value will be interpreted by the irradiance model.
%
%   The 25,122 km difference between them (~1.6%) is small relative to
%   constellation scales. It does not affect shading calculations.
%
% Authors: Planetary Sunshade Foundation
% Depends on: L1_Stability_Region_Data.xlsx (Glasgow / PSF partnership)

clear; clc;
fprintf('=== Glasgow Equilibrium Region Preprocessor ===\n\n');

% -------------------------------------------------------------------------
% USER SETTING — edit this path before running
% -------------------------------------------------------------------------
glasgow_file = '/Users/morgangoodwin/Desktop/PSF/MatLab - Planetary Sunshade Foundation/matlab/sunshade models/psf/excel/psf model/L1_Stability_Region_Data.xlsx';
output_file  = fullfile(fileparts(mfilename('fullpath')), 'equilibrium_envelope.mat');

% -------------------------------------------------------------------------
% CONSTANTS
% -------------------------------------------------------------------------
AU_km         = 149597870.7;       % 1 AU in km (nominal, IAU standard)
scale_factor  = 1.08;              % ellipse expansion: captures boundary
                                   % bulges of 5-15% beyond inscribed fit
safety_margin = 0.95;              % acceptance threshold for placement:
                                   % (Y/a)^2 + (Z/b)^2 <= safety_margin

% L1 (Glasgow): Hill sphere approximation (mu/3)^(1/3)
% Used to convert Glasgow AU coordinates — keeps envelope self-consistent
% with the Glasgow source data. Produces X values in [858, 866535] km.
mass_ratio    = 1 / 333000;
L1_frac       = (mass_ratio / 3)^(1/3);
L1_glasgow_km = L1_frac * AU_km;          % ~1,496,478 km
L1_x_AU       = 1 - L1_frac;             % synodic X of L1 per Glasgow

% L1 (SPICE): from PSF irradiance model (distance_PL1___Sr___Sr.m)
% Not used for coordinate conversion — stored for reference and for
% computing X_optimal_km, which the irradiance model will interpret.
L1_spice_km   = 1521600;                  % ~1,521,600 km

fprintf('Constants:\n');
fprintf('  AU                      : %.3f km\n', AU_km);
fprintf('  L1 Glasgow (Hill sphere): %.0f km\n', L1_glasgow_km);
fprintf('  L1 SPICE   (irrad model): %.0f km\n', L1_spice_km);
fprintf('  L1 difference           : %.0f km  (%.1f%%)\n', ...
        L1_spice_km - L1_glasgow_km, ...
        100 * (L1_spice_km - L1_glasgow_km) / L1_glasgow_km);
fprintf('  L1 synodic X (Glasgow)  : %.10f AU\n', L1_x_AU);
fprintf('  Ellipse scale factor    : %.2f\n', scale_factor);
fprintf('  Safety margin           : %.2f\n\n', safety_margin);

% -------------------------------------------------------------------------
% LOAD RAW GLASGOW DATA
% -------------------------------------------------------------------------
fprintf('Loading: %s\n', glasgow_file);

opts           = detectImportOptions(glasgow_file);
opts.DataRange = 'A2';
raw            = readmatrix(glasgow_file, opts);

X_AU_raw   = raw(:, 1);
Y_AU_raw   = raw(:, 2);
Z_AU_raw   = raw(:, 3);
Theta_raw  = raw(:, 4);

% Strip NaN rows (trailing blank rows from Excel range)
valid      = ~isnan(X_AU_raw) & ~isnan(Y_AU_raw) & ~isnan(Z_AU_raw);
X_AU_raw   = X_AU_raw(valid);
Y_AU_raw   = Y_AU_raw(valid);
Z_AU_raw   = Z_AU_raw(valid);
Theta_raw  = Theta_raw(valid);

fprintf('  Raw points loaded       : %d\n', sum(valid));
fprintf('  Glasgow X AU range      : [%.10f, %.10f]\n', ...
        min(X_AU_raw), max(X_AU_raw));

% -------------------------------------------------------------------------
% CONVERT TO L1-CENTRED KM
% -------------------------------------------------------------------------
% Glasgow X is synodic distance from the Sun (~0.984-0.990 AU).
% Sunward offset from L1 = (L1_x_AU - spacecraft_X_AU) * AU_km
% Positive X = sunward of L1, the only stable side for sail craft.
% Y and Z are already centred on zero and scale directly by AU_km.

X_km = (L1_x_AU - X_AU_raw) .* AU_km;
Y_km =  Y_AU_raw             .* AU_km;
Z_km =  Z_AU_raw             .* AU_km;

fprintf('\nConverted to L1-centred km:\n');
fprintf('  X range : [%9.1f, %9.1f] km  (sunward of L1)\n', min(X_km), max(X_km));
fprintf('  Y range : [%9.1f, %9.1f] km\n', min(Y_km), max(Y_km));
fprintf('  Z range : [%9.1f, %9.1f] km\n', min(Z_km), max(Z_km));

if min(X_km) < 0
    n_earthward = sum(X_km < 0);
    warning(['%d points have negative X (Earthward of L1). ' ...
             'This may indicate a mismatch between the Glasgow L1 ' ...
             'reference and the value used here. Review constants.'], ...
             n_earthward);
end

% -------------------------------------------------------------------------
% FIT ELLIPSE ENVELOPE  a(X), b(X), theta(X)
% -------------------------------------------------------------------------
% Bin data into X-slabs. In each slab:
%   a = max(|Y|)  — tangential semi-axis of the inscribed ellipse
%   b = max(|Z|)  — out-of-plane semi-axis of the inscribed ellipse
%   theta = median sail pitch angle of boundary-ish points (ellipse val > 0.7)
%           Falls back to median of all bin points if boundary is sparse.
%
% scale_factor expands the inscribed ellipse to capture the ~5-15% boundary
% bulge seen in the Glasgow data (largest at ~60/120 degree diagonals).
% theta is not scaled — it is a physical sail angle.

n_bins   = 80;
X_edges  = linspace(min(X_km), max(X_km), n_bins + 1);
X_mid    = 0.5 * (X_edges(1:end-1) + X_edges(2:end));

a_raw     = nan(n_bins, 1);
b_raw     = nan(n_bins, 1);
theta_raw = nan(n_bins, 1);
n_pts     = zeros(n_bins, 1);

for i = 1:n_bins
    in_bin   = X_km >= X_edges(i) & X_km < X_edges(i+1);
    n_pts(i) = sum(in_bin);
    if n_pts(i) >= 5
        a_raw(i) = max(abs(Y_km(in_bin)));
        b_raw(i) = max(abs(Z_km(in_bin)));

        % Theta: median over boundary-ish points in this bin
        a_approx = a_raw(i);
        b_approx = b_raw(i);
        Y_bin    = Y_km(in_bin);
        Z_bin    = Z_km(in_bin);
        T_bin    = Theta_raw(in_bin);
        ev       = (Y_bin / a_approx).^2 + (Z_bin / b_approx).^2;
        boundary = ev > 0.7;
        if sum(boundary) >= 3
            theta_raw(i) = median(T_bin(boundary));
        else
            theta_raw(i) = median(T_bin);
        end
    end
end

% Keep well-populated bins only
good      = ~isnan(a_raw) & ~isnan(b_raw) & n_pts >= 5;
X_env     = X_mid(good)';
a_env     = a_raw(good) .* scale_factor;
b_env     = b_raw(good) .* scale_factor;
theta_env = theta_raw(good);

% Light smoothing (3-pt moving average) to remove sampling noise
a_env     = movmean(a_env, 3);
b_env     = movmean(b_env, 3);
theta_env = movmean(theta_env, 3);

fprintf('\nEnvelope fit:\n');
fprintf('  Bins fitted             : %d / %d\n', sum(good), n_bins);
fprintf('  X coverage              : [%.0f, %.0f] km\n', min(X_env), max(X_env));
fprintf('  Peak a (Y semi-axis)    : %.0f km at X = %.0f km\n', ...
        max(a_env), X_env(a_env == max(a_env)));
fprintf('  Peak b (Z semi-axis)    : %.0f km at X = %.0f km\n', ...
        max(b_env), X_env(b_env == max(b_env)));

% -------------------------------------------------------------------------
% OPTIMAL ZONE
% -------------------------------------------------------------------------
% McInnes / Sanchez optimal distance: ~2.36M km from Earth.
% X_optimal_km uses L1_spice_km because the irradiance model (which will
% consume this value) uses the SPICE L1. The envelope geometry itself
% uses L1_glasgow_km — these are different things.
X_optimal_km      = 2360000 - L1_spice_km;
a_at_optimal      = interp1(X_env, a_env,     X_optimal_km, 'linear', 'extrap');
b_at_optimal      = interp1(X_env, b_env,     X_optimal_km, 'linear', 'extrap');
theta_at_optimal  = interp1(X_env, theta_env, X_optimal_km, 'linear', 'extrap');

fprintf('\nOptimal zone (2.36M km from Earth):\n');
fprintf('  X from L1 (SPICE-based) : %.0f km\n',        X_optimal_km);
fprintf('  Envelope a (Y)          : %.0f km\n',        a_at_optimal);
fprintf('  Envelope b (Z)          : %.0f km\n',        b_at_optimal);
fprintf('  Theta at optimal        : %.4f rad (%.2f deg)\n', ...
        theta_at_optimal, rad2deg(theta_at_optimal));

% -------------------------------------------------------------------------
% PACKAGE AND SAVE
% -------------------------------------------------------------------------
envelope.X_km             = X_env;
envelope.a_km             = a_env;
envelope.b_km             = b_env;
envelope.theta_rad        = theta_env;

envelope.scale_factor     = scale_factor;
envelope.safety_margin    = safety_margin;

envelope.L1_glasgow_km    = L1_glasgow_km;   % used for coordinate conversion
envelope.L1_spice_km      = L1_spice_km;     % used by irradiance model
envelope.AU_km            = AU_km;

envelope.X_optimal_km     = X_optimal_km;
envelope.a_at_optimal     = a_at_optimal;
envelope.b_at_optimal     = b_at_optimal;
envelope.theta_at_optimal = theta_at_optimal;

envelope.source_file      = glasgow_file;
envelope.date_created     = datestr(now, 'yyyy-mm-dd HH:MM:SS');

save(output_file, 'envelope');
fprintf('\nSaved: %s\n', output_file);

% -------------------------------------------------------------------------
% SANITY PLOTS
% -------------------------------------------------------------------------
figure('Name', 'Equilibrium Envelope Verification', 'NumberTitle', 'off', ...
       'Position', [100 100 1100 420]);

% Plot 1: a(X) and b(X) across full X range
subplot(1, 3, 1);
hold on;
plot(X_env / 1e3, a_env / 1e3, 'b-',  'LineWidth', 1.5, 'DisplayName', 'a  (Y semi-axis)');
plot(X_env / 1e3, b_env / 1e3, 'g-',  'LineWidth', 1.5, 'DisplayName', 'b  (Z semi-axis)');
xline(X_optimal_km / 1e3, 'r--', 'LineWidth', 1, 'DisplayName', 'Optimal zone');
xlabel('X from L1 (thousand km)');
ylabel('Semi-axis (thousand km)');
title('Envelope semi-axes vs X depth');
legend('Location', 'northwest');
grid on; box on;

% Plot 2: cross-section ellipse at optimal zone
subplot(1, 3, 2);
theta_plot = linspace(0, 2*pi, 300);
a_opt = a_at_optimal;
b_opt = b_at_optimal;
fill(a_opt * cos(theta_plot) / 1e3, ...
     b_opt * sin(theta_plot) / 1e3, ...
     [0.85 0.92 1.0], 'EdgeColor', [0.2 0.5 0.8], 'LineWidth', 1.5);
hold on;
r_safe = sqrt(safety_margin);
plot(a_opt * r_safe * cos(theta_plot) / 1e3, ...
     b_opt * r_safe * sin(theta_plot) / 1e3, ...
     'b--', 'LineWidth', 1, 'DisplayName', 'Safety margin');
xlabel('Y (thousand km)');
ylabel('Z (thousand km)');
title(sprintf('Cross-section at optimal zone\n(X = %.0f km from L1)', X_optimal_km));
legend({'Envelope', 'Safety margin (0.95)'}, 'Location', 'northeast');
axis equal; grid on; box on;

% Plot 3: theta(X) across full X range
subplot(1, 3, 3);
plot(X_env / 1e3, rad2deg(theta_env), 'm-', 'LineWidth', 1.5);
hold on;
xline(X_optimal_km / 1e3, 'r--', 'LineWidth', 1);
xlabel('X from L1 (thousand km)');
ylabel('Theta (degrees)');
title('Sail pitch angle \theta vs X depth');
grid on; box on;
