%% Evolution of L1 point during a timeframe

close all
clear all
clc

tic  % Start stopwatch timer

%% Pre-Set Constants

radius_Earth = 6371; % [km]
radius_Sun = 696340; % [km]
mass_Earth = 5.9722e24; % [kg]
mass_Sun = 1.988400e30; % [kg]
omega_Earth = 2*pi/(365.25*86400); % [rad/s]
G = 6.67430e-11; % [m^3/kg/s^2]
solar_luminosity = 3.846e26; % [W]
speed_light = 2.99792458e8; % [m/s]
sigma_SB = 5.670374419e-8; % W/m^2/K^4
T_sun   = 5772;             % K
R_sun   = 6.957e8;          % m

%% SPICE Initialization

path_to_MuSCAT_Supporting_Files = '../MuSCAT_Supporting_Files/';
path_to_MuSCAT_v2 = '../MuSCAT_Matlab_v2/';

addpath(genpath(path_to_MuSCAT_Supporting_Files));
addpath(genpath(path_to_MuSCAT_v2));
cspice_furnsh([path_to_MuSCAT_Supporting_Files,'SPICE/de440s.bsp']);
cspice_furnsh([path_to_MuSCAT_Supporting_Files,'SPICE/naif0012.tls']);
cspice_furnsh([path_to_MuSCAT_Supporting_Files,'SPICE/pck00011.tpc']);

%% Initialize Simulation

time_utc = '2025-01-03T12:00:00';
num_rays = 1e6;
num_years = 1; % Defines the number of years

all_rays_array = zeros(num_rays,14);

%% Earth Position and Orientation at reference time

et_t = cspice_str2et(time_utc);
state = cspice_spkezr('EARTH', et_t, 'J2000', 'NONE', 'SUN');
Earth_pos = state(1:3)'; % [km]
distance_Sun_Earth = norm(Earth_pos); % [km]

%% L1 Evolution Setup

% Define time range: 1 year, daily steps
% start_utc = '2025 JAN 01 00:00:00';
% stop_utc  = '2026 JAN 01 00:00:00';
% step_days = 1;
% et_start = cspice_str2et(start_utc);
% et_stop  = cspice_str2et(stop_utc);
% et_vec   = et_start : step_days*86400 : et_stop;
%

%% Initialize Time
t0_date = '01-JAN-2025 00:00:00'; % start date
t0 = 0;
time_step = 1*24*60*60;            % 1 day [s]
time_final = num_years*365*24*60*60;         % x years [s]

% Time Array
time_array = (t0:time_step:time_final)';

% Convert to SPICE Ephemeris Time
t0_time_date = cspice_str2et(t0_date);       % [sec since J2000]
time_date_array = t0_time_date + time_array; % [sec since J2000] for each day

% Datetime array for plotting
Matlab_date_array = datetime(2025,1,1) + days(time_array/86400);

% Preallocate arrays
N = length(time_date_array);
Earth_pos_J2000 = zeros(N,3);
omega_Earth_J2000 = zeros(N,1);
L1_pos_J2000 = zeros(N,3);
distance_Sun_Earth = zeros(N,1);
distance__L1 = zeros(N,1);
distance_Earth_L1 = zeros(N,1);
solar_constant_Earth = zeros(N,1); % [W/m^2]

%% Compute Evolution of SEL1
for i = 1:N
    et = time_date_array(i);

    % Earth heliocentric position (J2000 frame)
    state = cspice_spkezr('EARTH', et, 'J2000', 'NONE', 'SUN');
    Earth_pos = state(1:3)'; % [km]
    Earth_pos_J2000(i,:) = Earth_pos;

    % Earth Omega
    Earth_velocity = state(4:6)'; % [km/sec]

    Sun_pos_normalized = func_normalize_vec(-Earth_pos);
    Earth_velocity_normalized = func_normalize_vec(Earth_velocity);
    Sun_pos_Earth_velocity_angle = rad2deg(func_angle_between_vectors(Sun_pos_normalized', Earth_velocity_normalized')); % [rad]

    Earth_velocity_normal_plane = Earth_velocity*(sind(Sun_pos_Earth_velocity_angle));

    omega_Earth = norm(Earth_velocity_normal_plane)/norm(Earth_pos); % [rad/s]
    omega_Earth_J2000(i,1) = omega_Earth; % [rad/s]

    % Sun–Earth distance
    r_SE_km = norm(Earth_pos);
    r_SE_m  = r_SE_km * 1e3; % [m]

    % Solar constant on Earth
    solar_constant_Earth(i) = sigma_SB * T_sun^4 * (R_sun / r_SE_m)^2; % [W/m^2]

    % Solve for Sun–Earth L1
    syms x
    eqn = (G*mass_Sun/x^2) - (G*mass_Earth/(r_SE_m - x)^2) - (omega_Earth^2)*x == 0;
    S = double(solve(eqn, x));

    % Select valid positive solution
    x_L1_m = min(S(S > 0 & S < r_SE_m));
    x_L1_km = x_L1_m / 1e3;

    % Store distances
    distance_Sun_Earth(i) = r_SE_km;                     % Sun–Earth
    distance_Sun_L1(i)    = x_L1_km;                     % Sun–L1
    distance_Earth_L1(i)  = r_SE_km - x_L1_km;           % Earth–L1

    % L1 position vector (along Sun–Earth line)
    L1_pos_J2000(i,:) = (x_L1_km / r_SE_km) * Earth_pos;
end


%% Plot Results

% 3D Plot: Evolution of Sun–Earth L1 and Earth Orbit

sun_size_increase_factor = 40;
earth_size_increase_factor = 6000;

plot_name = figure();
clc
set(plot_name,'units','normalized','outerposition',[0 0 1 1])
hold on; grid on; axis equal;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Evolution of Sun–Earth L1 Point (J2000 Frame)');

% Plot trajectories
plot3(L1_pos_J2000(:,1), L1_pos_J2000(:,2), L1_pos_J2000(:,3), 'r', 'LineWidth', 1.5);
plot3(Earth_pos_J2000(:,1), Earth_pos_J2000(:,2), Earth_pos_J2000(:,3), 'b', 'LineWidth', 1.2);

% Textured Sun
sunTex = imread([path_to_MuSCAT_Supporting_Files,'SB_data/Sun/sun_texture.jpg']);
[sunTexRows, sunTexCols, ~] = size(sunTex);
if sunTexRows ~= sunTexCols
    sunTex = imresize(sunTex, [256 256]);
end

[sx, sy, sz] = sphere(100);
r_sun = 6.957e5 * sun_size_increase_factor; % scaled Sun radius [km]
hSun = surf(r_sun*sx, r_sun*sy, r_sun*sz, ...
    'FaceColor','texturemap', ...
    'CData', sunTex, ...
    'EdgeColor','none', ...
    'FaceLighting','none');

light('Position',[1 0 0],'Style','infinite');
material dull;

legend('Sun–Earth L1','Earth Orbit','Sun','Location','best');
view(45,20);

set(gca,'FontSize',14, 'FontName','Times New Roman')

saveas(plot_name, 'Sun_Earth_L1_Orbit.png');

%% Plot Distance Evolution: Earth–L1

plot_name = figure();
clc
set(plot_name,'units','normalized','outerposition',[0 0 1 1])

yyaxis left
p1 = plot(Matlab_date_array, distance_Earth_L1, '-', 'LineWidth', 2);
xlabel('Date');
ylabel('Earth–L_1 Distance [km]');

ylim_diff = max(distance_Earth_L1) - min(distance_Earth_L1);
ylim([min(distance_Earth_L1)-ylim_diff/10, max(distance_Earth_L1)+ylim_diff/10])
grid on;

yyaxis right
dist_frac = distance_Earth_L1 ./ min(distance_Earth_L1);
p2 = plot(Matlab_date_array, dist_frac, '--', 'LineWidth', 3);
ylabel('(Earth–L_1 Distance) / Min(Earth–L_1 Distance) ');
% legend([p1, p2], {'Earth–L_1 Distance', '(Earth–L_1) / Min(Earth–L_1)'}, 'Location', 'best');

% ylim([0.5 3]);
ylim_diff = max(dist_frac) - min(dist_frac);
ylim([min(dist_frac)-ylim_diff/10, max(dist_frac)+ylim_diff/10])


title('Evolution of Earth–L_1 Distance over the Year');

set(gca,'FontSize',14, 'FontName','Times New Roman')

saveas(plot_name, 'Earth_L1_Distance.png');

%% Plot Distance Evolution: Sun–L1

plot_name = figure();
clc
set(plot_name,'units','normalized','outerposition',[0 0 1 1])

yyaxis left
p3 = plot(Matlab_date_array, distance_Sun_L1, '-', 'LineWidth', 2); hold on;
% p4 = plot(Matlab_date_array, distance_Sun_Earth, '-k', 'LineWidth', 2);
xlabel('Date');
ylabel('Sun–L_1 Distance [km]');
grid on;

% Combine the datasets on left-axis
% all_data_left = [distance_Sun_L1(:); distance_Sun_Earth(:)];
ylim_diff = max(distance_Sun_L1) - min(distance_Sun_L1);
ylim([min(distance_Sun_L1) - ylim_diff/10, max(distance_Sun_L1) + ylim_diff/10]);

yyaxis right
dist_frac = distance_Sun_L1 ./ min(distance_Sun_L1);
p5 = plot(Matlab_date_array, dist_frac, '--', 'LineWidth', 3);
ylabel('(Sun–L_1 Distance) / Min(Sun–L_1 Distance)');

ylim_diff = max(dist_frac) - min(dist_frac);
ylim([min(dist_frac) - ylim_diff/10, max(dist_frac) + ylim_diff/10]);

% legend([p3, p4, p5], ...
%     {'Sun–L_1 Distance', 'Sun–Earth Distance', '(Sun–L_1) / Min(Sun–L_1)'}, ...
%     'Location', 'best');

title('Evolution of Sun–L_1 Distance over the Year');

set(gca, 'FontSize', 14, 'FontName', 'Times New Roman');

saveas(plot_name, 'Sun_L1_Distance.png');

%% Plot Distance Evolution: Sun–Earth

plot_name = figure();
clc
set(plot_name,'units','normalized','outerposition',[0 0 1 1])

yyaxis left
p3 = plot(Matlab_date_array, distance_Sun_Earth, '-', 'LineWidth', 2); hold on;
xlabel('Date');
ylabel('Sun–Earth Distance [km]');
grid on;

ylim_diff = max(distance_Sun_Earth) - min(distance_Sun_Earth);
ylim([min(distance_Sun_Earth) - ylim_diff/10, max(distance_Sun_Earth) + ylim_diff/10]);

yyaxis right
dist_frac = distance_Sun_Earth ./ min(distance_Sun_Earth);
p5 = plot(Matlab_date_array, dist_frac, '--', 'LineWidth', 3);
ylabel('(Sun–Earth Distance) / Min(Sun–Earth Distance)');

ylim_diff = max(dist_frac) - min(dist_frac);
ylim([min(dist_frac) - ylim_diff/10, max(dist_frac) + ylim_diff/10]);

% legend([p3, p4, p5], ...
%     {'Sun–L_1 Distance', 'Sun–Earth Distance', '(Sun–L_1) / Min(Sun–L_1)'}, ...
%     'Location', 'best');

title('Evolution of Sun–Earth Distance over the Year');

set(gca, 'FontSize', 14, 'FontName', 'Times New Roman');

saveas(plot_name, 'Sun_Earth_Distance.png');

%% Plot Solar Constant on Earth

plot_name = figure();
clc
set(plot_name,'units','normalized','outerposition',[0 0 1 1])

yyaxis left
p6 = plot(Matlab_date_array, solar_constant_Earth, '-', 'LineWidth', 2);
xlabel('Date');
ylabel('Solar Constant on Earth [W/m²]');
grid on;

ylim_diff = max(solar_constant_Earth) - min(solar_constant_Earth);
ylim([min(solar_constant_Earth) - ylim_diff/10, ...
    max(solar_constant_Earth) + ylim_diff/10]);


yyaxis right
dist_frac = solar_constant_Earth ./ min(solar_constant_Earth);
p5 = plot(Matlab_date_array, dist_frac, '--', 'LineWidth', 3);
ylabel('(Solar Constant) / Min(Solar Constant)');

ylim_diff = max(dist_frac) - min(dist_frac);
ylim([min(dist_frac) - ylim_diff/10, max(dist_frac) + ylim_diff/10]);

title('Evolution of the Solar Constant on Earth over the Year');

set(gca, 'FontSize', 14, 'FontName', 'Times New Roman');

saveas(plot_name, 'Solar_Constant_Earth.png');


%% Plot Earth Omega

plot_name = figure();
clc
set(plot_name,'units','normalized','outerposition',[0 0 1 1])

yyaxis left
p6 = plot(Matlab_date_array, omega_Earth_J2000, '-', 'LineWidth', 2);
xlabel('Date');
ylabel('Omega Earth [rad/sec]');
grid on;

ylim_diff = max(omega_Earth_J2000) - min(omega_Earth_J2000);
ylim([min(omega_Earth_J2000) - ylim_diff/10, ...
    max(omega_Earth_J2000) + ylim_diff/10]);


yyaxis right
dist_frac = omega_Earth_J2000 ./ min(omega_Earth_J2000);
p5 = plot(Matlab_date_array, dist_frac, '--', 'LineWidth', 3);
ylabel('(Solar Constant) / Min(Solar Constant)');

ylim_diff = max(dist_frac) - min(dist_frac);
ylim([min(dist_frac) - ylim_diff/10, max(dist_frac) + ylim_diff/10]);

title('Evolution of the Omega Earth over the Year');

set(gca, 'FontSize', 14, 'FontName', 'Times New Roman');

saveas(plot_name, 'Omega_Earth.png');



toc % Stop stopwatch timer

