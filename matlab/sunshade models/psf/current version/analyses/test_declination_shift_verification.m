%% test_declination_shift_verification.m
%
%  Verifies the solar declination coordinate correction by printing the
%  192-cell geographic IF profile for a set of key dates to the terminal.
%
%  The test can run in two modes:
%
%    REAL MODE   — load a pre-saved 192x1xN_days results array from a .mat
%                  file produced by a preconfigured run.  Save the raw
%                  (pre-shift) results inside analysis_General with:
%                      save(results_save_path, 'results')
%                  just after the analysis_Iterate_And_Switch call.
%
%    SYNTHETIC   — if no .mat file is found, a Gaussian-shaped IF profile
%                  is generated (equatorial constellation, no Z-offset).
%                  This is sufficient to verify the index arithmetic and
%                  the shape of the seasonal sweep.
%
%  The printed table shows, for each test date:
%    - Solar declination delta (deg) from Spencer (1971)
%    - Integer cell shift applied
%    - The full 192-row IF column, labelled with geographic latitude
%    - The row index of the minimum IF (expected to track delta)
%
%  Run this script directly from the MATLAB command window.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  USER SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Path to a pre-saved pre-shift results .mat file (192 x 1 x N_days).
%  Leave as '' to use the synthetic profile.
results_mat_path = 'constellation_B2_5293crafts_2026-04-03.mat';

%  Test dates: {label, day-of-year}
test_dates = { ...
    'Jan 01  (winter solstice region)',   1; ...
    'Mar 21  (vernal equinox)',          80; ...
    'Jun 21  (summer solstice)',        172; ...
    'Sep 22  (autumnal equinox)',       265; ...
    'Dec 21  (winter solstice)',        355; ...
};

%  Grid constants (must match preconfiguration)
n_phi        = 192;
lat_step_deg = 0.942408376963357;   % degrees per cell, CESM f09 grid
n_buffer     = 27;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  LOAD OR SYNTHESISE THE PRE-SHIFT RESULTS COLUMN  (192 x 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(results_mat_path) && isfile(results_mat_path)

    fprintf('\nLoading pre-shift results from:\n  %s\n', results_mat_path);
    S = load(results_mat_path, 'results');

    %  Use a single representative day (day 80 ≈ equinox) so the
    %  heliocentric shadow peak sits at row ~96-97 before any shift.
    equinox_day = 80;
    if size(S.results, 3) >= equinox_day
        base_column = S.results(:, 1, equinox_day);   % 192 x 1
    else
        base_column = S.results(:, 1, end);
    end
    source_label = sprintf('real results, day %d of loaded file', equinox_day);

else

    fprintf('\nNo results .mat file specified or found — using synthetic profile.\n');
    fprintf('To use real results: set results_mat_path at the top of this script.\n\n');

    %  Synthetic Gaussian IF profile: shadow centred on heliocentric phi=0,
    %  peak reduction ~1%% (target irradiance factor ~0.99 at equator).
    %  Values near limb are ~0.9915, matching typical Architecture A output.
    phi_centers  = linspace(-90 + lat_step_deg/2, 90 - lat_step_deg/2, n_phi)';
    sigma_phi    = 35;    % degrees — controls shadow footprint width
    base_column  = 1 - 0.01 * exp(-phi_centers.^2 / (2 * sigma_phi^2));
    %  Add a small limb floor so edge values are ~0.9915, not 1.0000
    limb_floor   = 0.0001 * (1 - cosd(phi_centers));
    base_column  = base_column - limb_floor;
    source_label = 'synthetic (Gaussian equatorial profile)';

end

%  Replicate across all test days — for synthetic mode the same profile is
%  used for every day; for real mode see note above.
n_test  = size(test_dates, 1);
results = repmat(base_column, [1, 1, 365]);   % 192 x 1 x 365


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  BUILD EXTENDED ARRAY WITH BUFFER CELLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Buffer cells filled with the per-day edge value (same logic as
%  analysis_General___0___X.m).  In this test the profile is identical
%  every day, so the buffer values are constant across the time dimension.

buffer_south    = repmat(results(1,   :, :), [n_buffer, 1, 1]);
buffer_north    = repmat(results(end, :, :), [n_buffer, 1, 1]);
results_extended = [buffer_south; results; buffer_north];   % 246 x 1 x 365


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  LATITUDE AXIS  (geographic, degrees north)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Cell centres: the preconfiguration interval is [-90-half_step, 90+half_step]
%  with n_phi cells, so the first centre lands at exactly -90 and the last
%  at exactly +90.  linspace is used here to avoid the floating-point
%  rounding that causes the colon operator to produce 190 instead of 192.
lat_geo   = linspace(-90, 90, n_phi)';   % 192 x 1, degrees north


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  PRINT VERIFICATION TABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n');
fprintf('============================================================\n');
fprintf('  Declination shift verification\n');
fprintf('  Source: %s\n', source_label);
fprintf('============================================================\n\n');

%  Gather shifted columns for all test dates first
shifted_cols  = zeros(n_phi, n_test);
deltas        = zeros(1, n_test);
shifts        = zeros(1, n_test);

for t = 1:n_test
    doy            = test_dates{t, 2};
    delta_deg      = spencer_declination(doy, 365);
    shift          = round(delta_deg / lat_step_deg);
    row_s          = n_buffer + 1 - shift;

    assert(row_s >= 1 && row_s + n_phi - 1 <= size(results_extended, 1), ...
        'Buffer overrun for doy=%d (shift=%d). Increase n_buffer.', doy, shift);

    shifted_cols(:, t) = results_extended(row_s : row_s + n_phi - 1, 1, doy);
    deltas(t)          = delta_deg;
    shifts(t)          = shift;
end

%  Header
fprintf('Row  Lat(geo)  ');
for t = 1:n_test
    fprintf('| %-22s ', test_dates{t,1}(1:min(22,end)));
end
fprintf('\n');

fprintf('     (deg)     ');
for t = 1:n_test
    fprintf('|  d=%+.1f  shift=%+d  ', deltas(t), shifts(t));
end
fprintf('\n');

fprintf('%s\n', repmat('-', 1, 15 + n_test * 26));

%  Rows
for i = 1:n_phi
    fprintf('%03d  %+7.2f°  ', i, lat_geo(i));
    for t = 1:n_test
        fprintf('|      %.6f        ', shifted_cols(i, t));
    end
    fprintf('\n');
end

fprintf('\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  SUMMARY: WHERE IS THE MINIMUM IF IN EACH COLUMN?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('--- Shadow peak (minimum IF) ---\n');
fprintf('%-28s  %8s  %8s  %8s\n', 'Date', 'Min IF', 'Row', 'Lat (geo)');
fprintf('%s\n', repmat('-', 1, 58));
for t = 1:n_test
    [min_val, min_row] = min(shifted_cols(:, t));
    fprintf('%-28s  %.6f  %8d  %+7.2f°\n', ...
        test_dates{t,1}, min_val, min_row, lat_geo(min_row));
end

fprintf('\nExpected: shadow peak latitude should track solar declination delta.\n');
fprintf('At equinox (delta=0): peak near lat=0 (rows 96-97).\n');
fprintf('At Jun solstice (delta=+23.4 deg): peak near lat=+23.4 (rows ~71-72).\n');
fprintf('At Dec solstice (delta=-23.4 deg): peak near lat=-23.4 (rows ~121-122).\n\n');
