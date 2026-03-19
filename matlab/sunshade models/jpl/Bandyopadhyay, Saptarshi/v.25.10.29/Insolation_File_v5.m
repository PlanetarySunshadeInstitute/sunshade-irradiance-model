% Solar Incidence on Earth with and without Dust Cloud

close all,
clear all, clc

tic % Start stopwatch timer - measures time for code exection

%% Pre-Set Constants

radius_Earth = 6371; % [km]

radius_Sun = 696340; % [km]

mass_Earth = 5.9722e24; % [kg]

mass_Sun = 1.988400e30; % [kg]

omega_Earth = 2*pi/(365.25*86400); % [rad/sec]

G = 6.67430e-11; % [N⋅m2⋅kg−2]

solar_luminosity = 3.846e26; % [W]

speed_light = 2.99792458e8; % [m/s]

%% SPICE Initialization

path_to_MuSCAT_Supporting_Files = '../MuSCAT_Supporting_Files/';
path_to_MuSCAT_v2 = '../MuSCAT_Matlab_v2/';

addpath(genpath(path_to_MuSCAT_Supporting_Files));
addpath(genpath(path_to_MuSCAT_v2));
cspice_furnsh([path_to_MuSCAT_Supporting_Files,'SPICE/de440s.bsp']);
cspice_furnsh([path_to_MuSCAT_Supporting_Files,'SPICE/naif0012.tls']);
cspice_furnsh([path_to_MuSCAT_Supporting_Files,'SPICE/pck00011.tpc']);

%% Initialize Simulation

% time_utc = '2025-12-21T12:00:00'; % this_time_utc
time_utc = '2025-01-01T12:00:00'; % Perihelion
% time_utc = '2025-07-04T12:00:00'; % Aphelion

num_rays = 1e7; % this_num_rays

% save(['all_data_',num2str(num_rays),'_',time_utc,'.mat'])
% save(replace(['all_data_',num2str(num_rays),'_',time_utc,'.mat'],':','_'))

all_rays_array = zeros(num_rays,14); % [Boolean hits Dust Cloud, Boolean hits Earth, [x1 y1 z1], [x2 y2 z2], [x3 y3 z3], [x4 y4 z4] ]

%% Earth Position and Orientation

% Convert UTC to ET
et_t = cspice_str2et(time_utc);

% Earth position at the defined time
state = cspice_spkezr('EARTH', et_t, 'J2000', 'NONE', 'SUN');
Earth_pos = state(1:3)'; % [km]

distance_Sun_Earth = norm(Earth_pos); % [km]

% Earth Omega
Earth_velocity = state(4:6)'; % [km/sec]

Sun_pos_normalized = func_normalize_vec(-Earth_pos);
Earth_velocity_normalized = func_normalize_vec(Earth_velocity);
Sun_pos_Earth_velocity_angle = rad2deg(func_angle_between_vectors(Sun_pos_normalized', Earth_velocity_normalized')); % [rad]

Earth_velocity_normal_plane = Earth_velocity*(sind(Sun_pos_Earth_velocity_angle));

omega_Earth = norm(Earth_velocity_normal_plane)/norm(Earth_pos); % [rad/s]


% Solve for Sun-Earth L1 point

syms x

eqn = (G*mass_Sun/((x^2) * (1e6))) - (G*mass_Earth/( (distance_Sun_Earth - x)^2 * (1e6) )) - ((omega_Earth^2) * x * 1e3) == 0;

S = solve(eqn,x);

distance_Sun_L1 = double(S(1)); % [km]

distance_Earth_L1 = (distance_Sun_Earth - distance_Sun_L1); % [km]

ratio_L1_Earth = distance_Sun_L1/distance_Sun_Earth;

% Earth orientation at the defined time
rotation_matrix_Earth = cspice_pxform('IAU_EARTH','J2000',et_t);

solar_constant_Earth = solar_luminosity/(4*pi*(distance_Sun_Earth^2)* 1e6); % [W/m2]

%% Outer radius at Dust Cloud location

flag_architecture = 'DimSun (Architecture E) 0.9percent reduction';

switch flag_architecture

    case 'DimSun (Architecture E) 0.9percent reduction'
        
        x_location = distance_Sun_L1; % [km] Location of Obstruction

        temp_x = (distance_Sun_Earth * radius_Earth)/(radius_Sun - radius_Earth); % [km]
        radius_location = (radius_Earth / temp_x ) * (temp_x + distance_Sun_Earth - x_location); % [km]

        radius_obstruction = 1.8555e+03; % [km] Architecture E: From Dust_Cloud_Mass_Marks_Equations.m, line 549. Reduction in blue ray intensity by 13.2485% for 0.9% reduction in solar intensity
        factor_dimming_obstruction = 1 - 0.132485;


    case 'DimSun (Architecture E) 1.16percent reduction'

        x_location = distance_Sun_L1; % [km] Location of Obstruction

        temp_x = (distance_Sun_Earth * radius_Earth)/(radius_Sun - radius_Earth); % [km]
        radius_location = (radius_Earth / temp_x ) * (temp_x + distance_Sun_Earth - x_location); % [km]

        radius_obstruction = 2.1206e+03; % [km] Architecture E: From Dust_Cloud_Mass_Marks_Equations.m, line 549. Reduction in blue ray intensity by 12.5398% for 1.16% reduction in solar intensity
        factor_dimming_obstruction = 1 - 0.125398;


    case 'Architecture A (from Dean)'

        distance_Earth_L1_Dean = 2.36e6; % [km]
        x_location = distance_Sun_Earth - distance_Earth_L1_Dean; % [km]

        temp_x = (distance_Sun_Earth * radius_Earth)/(radius_Sun - radius_Earth); % [km]
        radius_location = (radius_Earth / temp_x ) * (temp_x + distance_Sun_Earth - x_location); % [km]

        radius_obstruction = 2000; % [km] -> Architecture A: From Dean. Reduction in blue ray intensity by 50% because 5000SC 20km-radius are inside this region. 
        factor_dimming_obstruction = 0.125; % = 0.25*(pi*20^2)*5000/(pi*2000^2)


    case 'Architecture A (from Jeff)'

        x_location = distance_Sun_L1; % [km] Location of Obstruction

        temp_x = (distance_Sun_Earth * radius_Earth)/(radius_Sun - radius_Earth); % [km]
        radius_location = (radius_Earth / temp_x ) * (temp_x + distance_Sun_Earth - x_location); % [km]

        radius_obstruction = 970; % [km] -> Architecture A: From Jeff. Reduction in blue ray intensity by 100%
        factor_dimming_obstruction = 0;

    otherwise
        error('Should not reach here!')

end




%% Rotation Matrix for Ray

% Aligning the Primary Vector
primary_vector = func_normalize_vec([distance_Sun_Earth 0 0]);
desired_primary_vector = func_normalize_vec(Earth_pos);

% Compute rotation matrix
if 1 - dot(primary_vector, desired_primary_vector) <= 1e-10
    Rotation_matrix_ray = eye(3);
    warning('Rotation_matrix_ray is Identity!')

elseif 1 + dot(primary_vector, desired_primary_vector) <= 1e-10

    third_vector = func_normalize_vec([primary_vector(2:3), primary_vector(1)]);
    v = func_normalize_vec(cross(primary_vector, third_vector));

    Rotation_matrix_ray = -eye(3) + 2*(v')*v;
    warning('Rotation_matrix_ray is Badly computed!')

else
    % Rodrigues' formula assuming two unit and opposite vectors
    v = cross(primary_vector, desired_primary_vector);
    Rotation_matrix_ray = eye(3) + skew(v) + skew(v)^2*(1/(1+dot(primary_vector, desired_primary_vector)));

    % This is Rodrigue's formula; it gives a proper rotation matrix
    % Rot_primary = eye(3) + skew(v) + skew(v)^2*(1-dot(primary_vector, desired_primary_vector))/(norm(v)^2);
end



%% Ray Tracing

flag_Sun_ray_generator = 'limb darkening';

for i=1:1:num_rays

    if mod(100*i, num_rays) == 0
        disp(i)
    end

    switch flag_Sun_ray_generator

        case 'bands'

            % Select random location on Sun for the start of ray
            x1 = 0; % [km]
            r1 = rand*radius_Sun; % [km]
            th1 = rand*2*pi; % [rad]
            y1 = r1*cos(th1); % [km]
            z1 = r1*sin(th1); % [km]

            % Select random location on Sun-Earth L1 for ray
            x2 = x_location; % [km]
            r2 = rand*radius_location; % [km]
            th2 = rand*2*pi; % [rad]
            y2 = r2*cos(th2); % [km]
            z2 = r2*sin(th2); % [km]

        case 'uniformly random'

            % Select random location on Sun for the start of ray
            x1 = 0; % [km]

            flag_rand_variables = 0;
            while flag_rand_variables == 0
                yr1 = 2*(rand-0.5);
                zr1 = 2*(rand-0.5);
                if norm([yr1 zr1]) < 1
                    flag_rand_variables = 1;
                end
            end

            y1 = radius_Sun*yr1; % [km]
            z1 = radius_Sun*zr1; % [km]

            % Select random location on Sun-Earth L1 for ray
            x2 = x_location; % [km]

            flag_rand_variables = 0;
            while flag_rand_variables == 0
                yr2 = 2*(rand-0.5);
                zr2 = 2*(rand-0.5);
                if norm([yr2 zr2]) < 1
                    flag_rand_variables = 1;
                end
            end

            y2 = radius_location*yr2; % [km]
            z2 = radius_location*zr2; % [km]

        case 'limb darkening'

            x1 = 0; % [km] start at Sun center in x

            flag_rand_variables = 0;
            while flag_rand_variables == 0
                % Pick a random point in unit circle
                yr1 = 2*(rand-0.5);
                zr1 = 2*(rand-0.5);
                r_norm = norm([yr1 zr1]);

                if r_norm < 1

                    % % Reject brighter rays
                    % u = 0.5; % limb darkening coefficient
                    % mu = sqrt(1 - r_norm^2);
                    % intensity = 1 - u*(1 - mu); % brightness factor

                    % Second method: https://en.wikipedia.org/wiki/Limb_darkening
                    phi = asind(r_norm); % [deg]
                    a0 = 0.3;
                    a1 = 0.93;
                    a2 = -0.23;
                    intensity = a0 + (a1 * cosd(phi)) + (a2 * cosd(phi) * cosd(phi));

                    % Rejection sampling - keep rays probability proportional to intensity
                    if rand <= intensity
                        flag_rand_variables = 1;
                    else
                        flag_rand_variables = 0; % reject brighter rays at the limb
                    end

                end
            end


            % Scale to true solar radius
            y1 = radius_Sun*yr1;
            z1 = radius_Sun*zr1;

            % Ray at L1 remains same as uniform
            x2 = x_location; % [km]
            flag_rand_variables = 0;
            while flag_rand_variables == 0
                yr2 = 2*(rand-0.5);
                zr2 = 2*(rand-0.5);
                if norm([yr2 zr2]) < 1
                    flag_rand_variables = 1;
                end
            end
            y2 = radius_location*yr2; % [km]
            z2 = radius_location*zr2; % [km]

            % Save ray intensity for later use
            all_rays_array(i,15) = intensity;

        otherwise

            error('Ray generator is not defined!')

    end

    % Does the ray pass through Dust Cloud
    if norm([y2 z2]) <= radius_obstruction
        % Yes
        all_rays_array(i,1) = 1;
    end

    % Does the ray hit Earth?
    x3 = distance_Sun_Earth; % [km]

    y3 = (y2 - y1)*(x3 - x1)/(x2 - x1) + y1; % [km]
    z3 = (z2 - z1)*(x3 - x1)/(x2 - x1) + z1; % [km]

    all_rays_array(i,3:11) = [x1 y1 z1 x2 y2 z2 x3 y3 z3];

    if (norm([y3 z3])/radius_Earth) <= 1
        % Yes
        all_rays_array(i,2) = 1;

        % % Intersection with Earth's Sphere (using Matlab solve)
        % syms t
        % eqn2 = norm(t*[x2 y2 z2] + (1-t)*[x3 y3 z3] - [distance_Sun_Earth 0 0]) == radius_Earth;
        % S2 = solve(eqn2,t);
        %
        % if double(S2(1)) > 0
        %     t = double(S2(1));
        % elseif double(S2(2)) > 0
        %     t = double(S2(2));
        % else
        %     error('t is wrong!')
        % end

        % Intersection with Earth's Sphere (using quadratic equation)

        A = (x2 - x3)^2 + (y2 - y3)^2 + (z2 - z3)^2;
        B = 2*(x2 - x3)*(x3 - distance_Sun_Earth) + 2*(y2 - y3)*y3 + 2*(z2 - z3)*z3;
        C = (x3 - distance_Sun_Earth)^2 + y3^2 + z3^2 - radius_Earth^2;

        t = (-B + sqrt(B^2 - 4*A*C))/(2*A);
        if t < 0
            t = (-B - sqrt(B^2 - 4*A*C))/(2*A);
        end

        Sphere_Intersect_point_ray = t*[x2 y2 z2] + (1-t)*[x3 y3 z3];

        all_rays_array(i,12:14) = Sphere_Intersect_point_ray;

        % % Earth Intersect Point in J2000 frame
        % Sphere_Intersect_point_J2000 = (Rotation_matrix_ray * Sphere_Intersect_point_ray')'; % [km]
        %
        % Sphere_Intersect_point_relative_Earth = Sphere_Intersect_point_J2000 - Earth_pos; % [km]
        %
        % [radius, lon, lat] = cspice_reclat(rotation_matrix_Earth' * Sphere_Intersect_point_relative_Earth'); % [radius, longitude [rad], latitude [rad] ]
        %
        % longitude = rad2deg(lon); % [deg]
        % latitude = rad2deg(lat); % [deg]
        %
        % all_rays_array(i,3) = latitude; % [deg]
        % all_rays_array(i,4) = longitude; % [deg]

    end


end

disp(['Rays hitting Earth = ',num2str(100*sum(all_rays_array(:,2) == 1)/num_rays),' % ']) % Approx 46%

disp(['Rays hitting Earth and Obstruction = ',num2str(100*sum( (all_rays_array(:,2) == 1) & (all_rays_array(:,1) == 1) )/num_rays),' % ']) % Approx 17%

%% Find Latitude and Longitude

all_rays_lat_long = zeros(num_rays,2); % [Latitude Longitude]

for i = 1:1:num_rays

    if all_rays_array(i,2) == 1

        Sphere_Intersect_point_ray = all_rays_array(i,12:14);

        % Earth Intersect Point in J2000 frame
        Sphere_Intersect_point_J2000 = (Rotation_matrix_ray * Sphere_Intersect_point_ray')'; % [km]

        Sphere_Intersect_point_relative_Earth = Sphere_Intersect_point_J2000 - Earth_pos; % [km]

        [radius, lon, lat] = cspice_reclat(rotation_matrix_Earth' * Sphere_Intersect_point_relative_Earth'); % [radius, longitude [rad], latitude [rad] ]

        longitude = rad2deg(lon); % [deg]
        latitude = rad2deg(lat); % [deg]

        all_rays_lat_long(i,1) = latitude; % [deg]
        all_rays_lat_long(i,2) = longitude; % [deg]

    end

end

toc % Stop stopwatch timer

% Save rays data

% save(['all_data_',num2str(num_rays),'_',time_utc,'.mat'])
% save(replace(['all_data_',num2str(num_rays),'_',time_utc,'.mat'],':','_'))

%% Plot World Map

plot_handle = figure(1);
clf
set(plot_handle, 'Name','Ray Tracing')
set(plot_handle,'Color',[1 1 1]);
set(plot_handle,'units','normalized','outerposition',[0 0 0.5 0.5])
set(plot_handle,'PaperPositionMode','auto');

hold on

% Plot latitude/longitude of rays that hit Earth without passing through the dust cloud (red dots)

idx_Earth = logical( (all_rays_array(:,2) == 1) & (all_rays_array(:,1) == 0) );
plot(all_rays_lat_long(idx_Earth,2), all_rays_lat_long(idx_Earth,1), '.r')

% Plot latitude/longitude of rays that hit Earth AND passed through the dust cloud (blue dots)

idx_Earth = logical( (all_rays_array(:,2) == 1) & (all_rays_array(:,1) == 1) );
plot(all_rays_lat_long(idx_Earth,2), all_rays_lat_long(idx_Earth,1), '.b')

% Load coastline data for map overlay

load coastlines.mat
plot(coastlon, coastlat, 'k','LineWidth',2)

% Format map visualization properties

grid on
axis equal

xlim([-180 180])
ylim([-90 90])

xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
title(['Date = ',time_utc])
set(gca,'FontSize',15, 'FontName','Times New Roman')

% Save image

% saveas(plot_handle,replace(['Rays_',num2str(num_rays),'_',time_utc,'.png'],':','_'))

%% Dust vs non-dust rays on Mercator map with 10x10° grid
fprintf('\n Dust vs non-dust rays on Mercator map with 10x10° grid \n');

idx_hit = (all_rays_array(:,2) == 1);
idx_dust = (all_rays_array(:,1) == 1);
idx_hit_dust = idx_hit & idx_dust;
idx_hit_nodust = idx_hit & ~idx_dust;
idx_hit_total = idx_hit;

% Extract coordinates
lat_hit = all_rays_lat_long(idx_hit,1);
lon_hit = all_rays_lat_long(idx_hit,2);
lat_dust = all_rays_lat_long(idx_hit_dust,1);
lon_dust = all_rays_lat_long(idx_hit_dust,2);
lat_nodust = all_rays_lat_long(idx_hit_nodust,1);
lon_nodust = all_rays_lat_long(idx_hit_nodust,2);
lat_total = all_rays_lat_long(idx_hit_total,1);
lon_total = all_rays_lat_long(idx_hit_total,2);

% Define 10x10° grid edges and centers
delta_lon = 10; % [deg]
delta_lat = 10; % [deg]
lon_edges = -180:delta_lon:180;
lat_edges = -90:delta_lat:90;
lon_centers = lon_edges(1:end-1) + diff(lon_edges)/2;
lat_centers = lat_edges(1:end-1) + diff(lat_edges)/2;

% Ray count
[counts_nodust, ~, ~] = histcounts2(lat_nodust, lon_nodust, lat_edges, lon_edges);
[counts_dust, ~, ~]  = histcounts2(lat_dust, lon_dust, lat_edges, lon_edges);
[counts_total, ~, ~]  = histcounts2(lat_total, lon_total, lat_edges, lon_edges);

% Load background map
I = imread('earth.png');

% Create figure
fig3 = figure('Name','Ray Tracing with Dust Cloud on Mercator Map','Color',[1 1 1]);
set(fig3,'units','normalized','outerposition',[0.4 0 0.6 0.65])
ax = gca; hold on

% Background map
if ~isempty(I)
    imagesc([-180 180], [-90 90], flipud(I));
    set(ax,'YDir','normal');
else
    xlim([-180 180]); ylim([-90 90]);
end

% Overlay heatmap of total rays per cell
h_counts = imagesc([-180 180], [-90 90], counts_nodust);
set(ax,'YDir','normal');
alpha(h_counts, 0.35);
colormap(parula)
colorbar
caxis([0 max(counts_nodust(:))]);
% title(sprintf('Solar Rays hitting Earth with a Dust Cloud (%s)', time_utc), 'FontWeight','bold')
title(['Date = ',time_utc,', Solar Rays hitting Earth with an Obstruction of Radius = ',num2str(radius_obstruction),' km at SEL_1'], 'FontWeight','bold')
xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
axis equal
xlim([-180 180]); ylim([-90 90]);

% Plot the rays
plot(lon_nodust, lat_nodust, '.r', 'MarkerSize', 4);
plot(lon_dust, lat_dust, '.b', 'MarkerSize', 4);

% Overlay 10° white grid
for lonV = lon_edges
    plot([lonV lonV], [-90 90], 'w', 'LineWidth', 0.5);
end
for latH = lat_edges
    plot([-180 180], [latH latH], 'w', 'LineWidth', 0.5);
end

% Annotate number of rays per cell (counts_total)

for r = 1:length(lat_centers)
    for c = 1:length(lon_centers)
        n = counts_nodust(r,c);
        if n > 0 % Zeros are skipped
            ndust = counts_dust(r,c);
            text(lon_centers(c), lat_centers(r), sprintf('%d', n), ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                'FontSize',7, 'FontWeight','bold', 'Color','w', ...
                'BackgroundColor',[0 0 0 0.4], 'Margin',1);
        end
    end
end

legend({'No Obstruction','Through Obstruction'}, 'Location','southoutside', 'Orientation','horizontal');
grid on
set(gca,'FontSize',14,'FontName','Times New Roman');

% Save figure
saveas(fig3, replace(['Mercator_Rays_',num2str(num_rays),'_Radius_',num2str(round(radius_obstruction)),'_Date_',time_utc,'.png'],':','_'));

% Print summary
fprintf('Total rays that hit Earth: %d\n', sum(idx_hit));
fprintf('  → Through Obstruction: %d\n', sum(idx_hit_dust));
fprintf('  → Without Obstruction: %d\n', sum(idx_hit_nodust));

%% Plot Intensity on Sphere

area_grid_cells_sphere = 0*counts_total;

for r = 1:length(lat_centers)
    for c = 1:length(lon_centers)

        this_lat_lower_edge = lat_edges(r);
        this_lat_upper_edge = lat_edges(r+1);

        area_grid_cells_sphere(r,c) = radius_Earth^2 * deg2rad(delta_lon) * (sind(this_lat_upper_edge) - sind(this_lat_lower_edge) ) * 1e6; % [m^2]

    end
end

% sum(sum(area_grid_cells)) = 4*pi*radius_Earth^2 = 5.1006e+14 km^2

total_energy_from_Sun_to_Earth = solar_constant_Earth * pi * radius_Earth^2 * 1e6; % [W]

energy_per_ray = total_energy_from_Sun_to_Earth/sum(idx_hit_total);  % [W]

intensity_total = energy_per_ray * counts_total ./area_grid_cells_sphere; % [W/m^2]

intensity_SRM = (energy_per_ray * counts_nodust ./area_grid_cells_sphere) + (energy_per_ray * factor_dimming_obstruction * counts_dust ./area_grid_cells_sphere); % [W/m^2]

intensity_reduction_SRM = intensity_total - intensity_SRM; % [W/m^2]

fraction_intensity_reduction_SRM = 100*intensity_reduction_SRM ./intensity_total; % [percentage]

% Create figure
fig3 = figure('Name','Solar Intensity without Obstruction','Color',[1 1 1]);
clc
set(fig3,'units','normalized','outerposition',[0.4 0 0.6 0.65])
ax = gca; hold on

% Background map
if ~isempty(I)
    imagesc([-180 180], [-90 90], flipud(I));
    set(ax,'YDir','normal');
else
    xlim([-180 180]); ylim([-90 90]);
end

% Overlay heatmap of total rays per cell
h_counts = imagesc([-180 180], [-90 90], intensity_total);
set(ax,'YDir','normal');
alpha(h_counts, 0.35);
colormap(parula)
a = colorbar;
a.Label.String = 'Solar Intensity [W/m^2]';
caxis([0 max(intensity_total(:))]);

% title(sprintf('Solar Rays hitting Earth with a Dust Cloud (%s)', time_utc), 'FontWeight','bold')
title(['Date = ',time_utc,', Solar Intensity on Earth, without any Obstruction'], 'FontWeight','bold')
xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
axis equal
xlim([-180 180]); ylim([-90 90]);

% Overlay 10° white grid
for lonV = lon_edges
    plot([lonV lonV], [-90 90], 'w', 'LineWidth', 0.5);
end
for latH = lat_edges
    plot([-180 180], [latH latH], 'w', 'LineWidth', 0.5);
end

% Annotate number of rays per cell (counts_total)

for r = 1:length(lat_centers)
    for c = 1:length(lon_centers)
        n = round(intensity_total(r,c));
        if n > 0 % Zeros are skipped
            ndust = counts_dust(r,c);
            text(lon_centers(c), lat_centers(r), sprintf('%d', n), ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                'FontSize',7, 'FontWeight','bold', 'Color','w', ...
                'BackgroundColor',[0 0 0 0.4], 'Margin',1);
        end
    end
end

set(gca,'FontSize',14, 'FontName','Times New Roman')

% Save figure
saveas(fig3, replace(['Intensity_no_Obstruction_Rays',num2str(num_rays),'_Radius_',num2str(round(radius_obstruction)),'_Date_',time_utc,'.png'],':','_'));


% Create figure
fig3 = figure('Name','Solar Intensity with SRM Obstruction','Color',[1 1 1]);
clc
set(fig3,'units','normalized','outerposition',[0 0 0.6 0.65])
ax = gca; hold on

% Background map
if ~isempty(I)
    imagesc([-180 180], [-90 90], flipud(I));
    set(ax,'YDir','normal');
else
    xlim([-180 180]); ylim([-90 90]);
end

% Overlay heatmap of total rays per cell
h_counts = imagesc([-180 180], [-90 90], intensity_SRM);
set(ax,'YDir','normal');
alpha(h_counts, 0.35);
colormap(parula)
a = colorbar;
a.Label.String = 'Solar Intensity [W/m^2]';
caxis([0 max(intensity_SRM(:))]);

% title(sprintf('Solar Rays hitting Earth with a Dust Cloud (%s)', time_utc), 'FontWeight','bold')
title(['Date = ',time_utc,', Solar Intensity on Earth with an Obstruction of Radius = ',num2str(radius_obstruction),' km at SEL_1'], 'FontWeight','bold')
xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
axis equal
xlim([-180 180]); ylim([-90 90]);

% Overlay 10° white grid
for lonV = lon_edges
    plot([lonV lonV], [-90 90], 'w', 'LineWidth', 0.5);
end
for latH = lat_edges
    plot([-180 180], [latH latH], 'w', 'LineWidth', 0.5);
end

% Annotate number of rays per cell (counts_total)

for r = 1:length(lat_centers)
    for c = 1:length(lon_centers)
        n = round(intensity_SRM(r,c));
        if n > 0 % Zeros are skipped
            ndust = counts_dust(r,c);
            text(lon_centers(c), lat_centers(r), sprintf('%d', n), ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                'FontSize',7, 'FontWeight','bold', 'Color','w', ...
                'BackgroundColor',[0 0 0 0.4], 'Margin',1);
        end
    end
end

set(gca,'FontSize',14, 'FontName','Times New Roman')

% Save figure
saveas(fig3, replace(['Intensity_with_SRM_Obstruction_Rays',num2str(num_rays),'_Radius_',num2str(round(radius_obstruction)),'_Date_',time_utc,'.png'],':','_'));


% Create figure
fig3 = figure('Name','Solar Intensity Reduction due to SRM Obstruction','Color',[1 1 1]);
clc
set(fig3,'units','normalized','outerposition',[0 0.3 0.6 0.65])
ax = gca; hold on

% Background map
if ~isempty(I)
    imagesc([-180 180], [-90 90], flipud(I));
    set(ax,'YDir','normal');
else
    xlim([-180 180]); ylim([-90 90]);
end

% Overlay heatmap of total rays per cell
h_counts = imagesc([-180 180], [-90 90], intensity_reduction_SRM);
set(ax,'YDir','normal');
alpha(h_counts, 0.35);
colormap(parula)
a = colorbar;
a.Label.String = 'Solar Intensity Reduction [W/m^2]';
caxis([0 max(intensity_reduction_SRM(:))]);

% title(sprintf('Solar Rays hitting Earth with a Dust Cloud (%s)', time_utc), 'FontWeight','bold')
title(['Date = ',time_utc,', Solar Intensity Reduction due to an Obstruction of Radius = ',num2str(radius_obstruction),' km at SEL_1'], 'FontWeight','bold')
xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
axis equal
xlim([-180 180]); ylim([-90 90]);

% Overlay 10° white grid
for lonV = lon_edges
    plot([lonV lonV], [-90 90], 'w', 'LineWidth', 0.5);
end
for latH = lat_edges
    plot([-180 180], [latH latH], 'w', 'LineWidth', 0.5);
end

% Annotate number of rays per cell (counts_total)

for r = 1:length(lat_centers)
    for c = 1:length(lon_centers)
        n = round(intensity_reduction_SRM(r,c));
        if n > 0 % Zeros are skipped
            ndust = counts_dust(r,c);
            text(lon_centers(c), lat_centers(r), sprintf('%d', n), ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                'FontSize',7, 'FontWeight','bold', 'Color','w', ...
                'BackgroundColor',[0 0 0 0.4], 'Margin',1);
        end
    end
end

set(gca,'FontSize',14, 'FontName','Times New Roman')

% Save figure
saveas(fig3, replace(['Intensity_Reduction_with_SRM_Obstruction_Rays',num2str(num_rays),'_Radius_',num2str(round(radius_obstruction)),'_Date_',time_utc,'.png'],':','_'));


% Create figure
fig3 = figure('Name','Fraction Solar Intensity Reduction due to SRM Obstruction','Color',[1 1 1]);
clc
set(fig3,'units','normalized','outerposition',[0.4 0.3 0.6 0.65])
ax = gca; hold on

% Background map
if ~isempty(I)
    imagesc([-180 180], [-90 90], flipud(I));
    set(ax,'YDir','normal');
else
    xlim([-180 180]); ylim([-90 90]);
end

% Overlay heatmap of total rays per cell
h_counts = imagesc([-180 180], [-90 90], fraction_intensity_reduction_SRM);
set(ax,'YDir','normal');
alpha(h_counts, 0.35);
colormap(parula)
a = colorbar;
a.Label.String = 'Fraction Solar Intensity Reduction [%]';
caxis([0 max(fraction_intensity_reduction_SRM(:))]);

% title(sprintf('Solar Rays hitting Earth with a Dust Cloud (%s)', time_utc), 'FontWeight','bold')
title(['Date = ',time_utc,', Fraction (%) Solar Intensity Reduction due to an Obstruction of Radius = ',num2str(radius_obstruction),' km at SEL_1'], 'FontWeight','bold')
xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
axis equal
xlim([-180 180]); ylim([-90 90]);

% Overlay 10° white grid
for lonV = lon_edges
    plot([lonV lonV], [-90 90], 'w', 'LineWidth', 0.5);
end
for latH = lat_edges
    plot([-180 180], [latH latH], 'w', 'LineWidth', 0.5);
end

% Annotate number of rays per cell (counts_total)

for r = 1:length(lat_centers)
    for c = 1:length(lon_centers)
        n = round(fraction_intensity_reduction_SRM(r,c),1);
        if n > 0 % Zeros are skipped
            ndust = counts_dust(r,c);
            text(lon_centers(c), lat_centers(r), sprintf('%s', num2str(n)), ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                'FontSize',7, 'FontWeight','bold', 'Color','w', ...
                'BackgroundColor',[0 0 0 0.4], 'Margin',1);
        end
    end
end

set(gca,'FontSize',14, 'FontName','Times New Roman')

% Save figure
saveas(fig3, replace(['Fraction_Intensity_Reduction_with_SRM_Obstruction_Rays',num2str(num_rays),'_Radius_',num2str(round(radius_obstruction)),'_Date_',time_utc,'.png'],':','_'));


%% Sun-Earth Ray Intersection

% Select center location on Sun for the start of ray
x1 = 0; % [km]
y1 = 0; % [km]
z1 = 0; % [km]

% Select random location on Sun-Earth L1 for ray
x2 = x_location; % [km]
y2 = 0; % [km]
z2 = 0; % [km]

% Does the ray hit Earth?
x3 = distance_Sun_Earth; % [km]

y3 = (y2 - y1)*(x3 - x1)/(x2 - x1) + y1; % [km]
z3 = (z2 - z1)*(x3 - x1)/(x2 - x1) + z1; % [km]

this_Sun_Earth_center_ray = [x1 y1 z1 x2 y2 z2 x3 y3 z3];

% Intersection with Earth's Sphere (using quadratic equation)
A = (x2 - x3)^2 + (y2 - y3)^2 + (z2 - z3)^2;
B = 2*(x2 - x3)*(x3 - distance_Sun_Earth) + 2*(y2 - y3)*y3 + 2*(z2 - z3)*z3;
C = (x3 - distance_Sun_Earth)^2 + y3^2 + z3^2 - radius_Earth^2;

t = (-B + sqrt(B^2 - 4*A*C))/(2*A);
if t < 0
    t = (-B - sqrt(B^2 - 4*A*C))/(2*A);
end

this_Sphere_Intersect_point_ray = t*[x2 y2 z2] + (1-t)*[x3 y3 z3];

% Earth Intersect Point in J2000 frame
this_Sphere_Intersect_point_J2000 = (Rotation_matrix_ray * this_Sphere_Intersect_point_ray')'; % [km]

this_Sphere_Intersect_point_relative_Earth = this_Sphere_Intersect_point_J2000 - Earth_pos; % [km]

[radius, lon, lat] = cspice_reclat(rotation_matrix_Earth' * this_Sphere_Intersect_point_relative_Earth'); % [radius, longitude [rad], latitude [rad] ]

Sun_longitude = rad2deg(lon); % [deg]
Sun_latitude = rad2deg(lat); % [deg]

%% Area Grid Cells on Projected Disc

area_grid_cells_projected_disc = 0*counts_total;

for r = 1:length(lat_centers)
    for c = 1:length(lon_centers)

        if counts_total(r,c) > 0

            % Inputs (all angles in DEGREES):
            %   phi1_deg, phi2_deg     - cell latitude bounds  (south, north)
            %   lambda1_deg, lambda2_deg - cell longitude bounds (west, east)
            %   phi0_deg, lambda0_deg  - subsolar latitude and longitude
            %   R                      - planetary radius (same units as output area)

            % Convert to radians
            d2r   = pi/180;
            phi1  = lat_edges(r)   * d2r;
            phi2  = lat_edges(r+1) * d2r;
            lam1  = lon_edges(c)   * d2r;
            lam2  = lon_edges(c+1) * d2r;
            phi0  = Sun_latitude   * d2r;
            lam0  = Sun_longitude  * d2r;

            R = radius_Earth; % [km]

            % Convenience deltas
            dlam = lam2 - lam1;

            % ---- Analytic integral over the cell ----
            % Integrand: mu * dA =
            %   (sinφ sinφ0 + cosφ cosφ0 cos(λ - λ0)) * R^2 cosφ dφ dλ
            %
            % Result:
            % Aproj = (R^2/2) * sinφ0 * dλ * (sin^2φ2 - sin^2φ1)
            %       + R^2 * cosφ0 * [sin(λ2 - λ0) - sin(λ1 - λ0)] ...
            %         * [ (φ2 - φ1)/2 + (sin2φ2 - sin2φ1)/4 ]

            term1 = 0.5 * R^2 * sin(phi0) * dlam * (sin(phi2)^2 - sin(phi1)^2);

            term2_lam = cos(phi0) * ( sin(lam2 - lam0) - sin(lam1 - lam0) );
            term2_phi = 0.5*(phi2 - phi1) + 0.25*(sin(2*phi2) - sin(2*phi1));

            Aproj = term1 + R^2 * term2_lam * term2_phi;

            if Aproj > 0
                area_grid_cells_projected_disc(r,c) = Aproj * 1e6; % [m^2]
            else
                area_grid_cells_projected_disc(r,c) = NaN;
            end


        else
            area_grid_cells_projected_disc(r,c) = NaN;

        end

    end
end

% Remove cells will very small area
max_grid_cell_area = max(max(area_grid_cells_projected_disc));
for r = 1:length(lat_centers)
    for c = 1:length(lon_centers)
        if area_grid_cells_projected_disc(r,c) < max_grid_cell_area/100
            area_grid_cells_projected_disc(r,c) = NaN;
        end
    end
end

%% Plot Intensity on Sun Projected Disc

intensity_total = energy_per_ray * counts_total ./area_grid_cells_projected_disc; % [W/m^2]

intensity_SRM = (energy_per_ray * counts_nodust ./area_grid_cells_projected_disc) + (energy_per_ray * factor_dimming_obstruction * counts_dust ./area_grid_cells_projected_disc); % [W/m^2]

intensity_reduction_SRM = intensity_total - intensity_SRM; % [W/m^2]

fraction_intensity_reduction_SRM = 100*intensity_reduction_SRM ./intensity_total; % [percentage]

% Create figure
fig3 = figure('Name','Solar Intensity on Sun Projected Disc without Obstruction','Color',[1 1 1]);
clc
set(fig3,'units','normalized','outerposition',[0.4 0 0.6 0.65])
ax = gca; hold on

% Background map
if ~isempty(I)
    imagesc([-180 180], [-90 90], flipud(I));
    set(ax,'YDir','normal');
else
    xlim([-180 180]); ylim([-90 90]);
end

% Overlay heatmap of total rays per cell
h_counts = imagesc([-180 180], [-90 90], intensity_total);
set(ax,'YDir','normal');
alpha(h_counts, 0.35);
colormap(parula)
a = colorbar;
a.Label.String = 'Solar Intensity on Projected Disc [W/m^2]';
caxis([0 max(intensity_total(:))]);

% title(sprintf('Solar Rays hitting Earth with a Dust Cloud (%s)', time_utc), 'FontWeight','bold')
title(['Date = ',time_utc,', Solar Intensity on Sun Projected Disc, without any Obstruction'], 'FontWeight','bold')
xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
axis equal
xlim([-180 180]); ylim([-90 90]);

% Overlay 10° white grid
for lonV = lon_edges
    plot([lonV lonV], [-90 90], 'w', 'LineWidth', 0.5);
end
for latH = lat_edges
    plot([-180 180], [latH latH], 'w', 'LineWidth', 0.5);
end

% Annotate number of rays per cell (counts_total)

for r = 1:length(lat_centers)
    for c = 1:length(lon_centers)
        n = round(intensity_total(r,c));
        if n > 0 % Zeros are skipped
            ndust = counts_dust(r,c);
            text(lon_centers(c), lat_centers(r), sprintf('%d', n), ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                'FontSize',7, 'FontWeight','bold', 'Color','w', ...
                'BackgroundColor',[0 0 0 0.4], 'Margin',1);
        end
    end
end

set(gca,'FontSize',14, 'FontName','Times New Roman')

% Save figure
saveas(fig3, replace(['Intensity_Projected_Disc_no_Obstruction_Rays',num2str(num_rays),'_Radius_',num2str(round(radius_obstruction)),'_Date_',time_utc,'.png'],':','_'));


% Create figure
fig3 = figure('Name','Solar Intensity on Sun Projected Disc with SRM Obstruction','Color',[1 1 1]);
clc
set(fig3,'units','normalized','outerposition',[0 0 0.6 0.65])
ax = gca; hold on

% Background map
if ~isempty(I)
    imagesc([-180 180], [-90 90], flipud(I));
    set(ax,'YDir','normal');
else
    xlim([-180 180]); ylim([-90 90]);
end

% Overlay heatmap of total rays per cell
h_counts = imagesc([-180 180], [-90 90], intensity_SRM);
set(ax,'YDir','normal');
alpha(h_counts, 0.35);
colormap(parula)
a = colorbar;
a.Label.String = 'Solar Intensity [W/m^2]';
caxis([0 max(intensity_SRM(:))]);

% title(sprintf('Solar Rays hitting Earth with a Dust Cloud (%s)', time_utc), 'FontWeight','bold')
title(['Date = ',time_utc,', Solar Intensity on Sun Projected Disc with an Obstruction of Radius = ',num2str(radius_obstruction),' km at SEL_1'], 'FontWeight','bold')
xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
axis equal
xlim([-180 180]); ylim([-90 90]);

% Overlay 10° white grid
for lonV = lon_edges
    plot([lonV lonV], [-90 90], 'w', 'LineWidth', 0.5);
end
for latH = lat_edges
    plot([-180 180], [latH latH], 'w', 'LineWidth', 0.5);
end

% Annotate number of rays per cell (counts_total)

for r = 1:length(lat_centers)
    for c = 1:length(lon_centers)
        n = round(intensity_SRM(r,c));
        if n > 0 % Zeros are skipped
            ndust = counts_dust(r,c);
            text(lon_centers(c), lat_centers(r), sprintf('%d', n), ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                'FontSize',7, 'FontWeight','bold', 'Color','w', ...
                'BackgroundColor',[0 0 0 0.4], 'Margin',1);
        end
    end
end

set(gca,'FontSize',14, 'FontName','Times New Roman')

% Save figure
saveas(fig3, replace(['Intensity_Projected_Disc_with_SRM_Obstruction_Rays',num2str(num_rays),'_Radius_',num2str(round(radius_obstruction)),'_Date_',time_utc,'.png'],':','_'));

% Create figure
fig3 = figure('Name','Solar Intensity Reduction on Sun Projected Disc due to SRM Obstruction','Color',[1 1 1]);
clc
set(fig3,'units','normalized','outerposition',[0 0.3 0.6 0.65])
ax = gca; hold on

% Background map
if ~isempty(I)
    imagesc([-180 180], [-90 90], flipud(I));
    set(ax,'YDir','normal');
else
    xlim([-180 180]); ylim([-90 90]);
end

% Overlay heatmap of total rays per cell
h_counts = imagesc([-180 180], [-90 90], intensity_reduction_SRM);
set(ax,'YDir','normal');
alpha(h_counts, 0.35);
colormap(parula)
a = colorbar;
a.Label.String = 'Solar Intensity Reduction [W/m^2]';
caxis([0 max(intensity_reduction_SRM(:))]);

% title(sprintf('Solar Rays hitting Earth with a Dust Cloud (%s)', time_utc), 'FontWeight','bold')
title(['Date = ',time_utc,', Solar Intensity Reduction on Sun Projected Disc due to an Obstruction of Radius = ',num2str(radius_obstruction),' km at SEL_1'], 'FontWeight','bold')
xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
axis equal
xlim([-180 180]); ylim([-90 90]);

% Overlay 10° white grid
for lonV = lon_edges
    plot([lonV lonV], [-90 90], 'w', 'LineWidth', 0.5);
end
for latH = lat_edges
    plot([-180 180], [latH latH], 'w', 'LineWidth', 0.5);
end

% Annotate number of rays per cell (counts_total)

for r = 1:length(lat_centers)
    for c = 1:length(lon_centers)
        n = round(intensity_reduction_SRM(r,c));
        if n > 0 % Zeros are skipped
            ndust = counts_dust(r,c);
            text(lon_centers(c), lat_centers(r), sprintf('%d', n), ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                'FontSize',7, 'FontWeight','bold', 'Color','w', ...
                'BackgroundColor',[0 0 0 0.4], 'Margin',1);
        end
    end
end

set(gca,'FontSize',14, 'FontName','Times New Roman')

% Save figure
saveas(fig3, replace(['Intensity_Reduction_Projected_Disc_with_SRM_Obstruction_Rays',num2str(num_rays),'_Radius_',num2str(round(radius_obstruction)),'_Date_',time_utc,'.png'],':','_'));


% Create figure
fig3 = figure('Name','Fraction Solar Intensity Reduction on Sun Projected Disc due to SRM Obstruction','Color',[1 1 1]);
clc
set(fig3,'units','normalized','outerposition',[0.4 0.3 0.6 0.65])
ax = gca; hold on

% Background map
if ~isempty(I)
    imagesc([-180 180], [-90 90], flipud(I));
    set(ax,'YDir','normal');
else
    xlim([-180 180]); ylim([-90 90]);
end

% Overlay heatmap of total rays per cell
h_counts = imagesc([-180 180], [-90 90], fraction_intensity_reduction_SRM);
set(ax,'YDir','normal');
alpha(h_counts, 0.35);
colormap(parula)
a = colorbar;
a.Label.String = 'Fraction Solar Intensity Reduction [%]';
caxis([0 max(fraction_intensity_reduction_SRM(:))]);

% title(sprintf('Solar Rays hitting Earth with a Dust Cloud (%s)', time_utc), 'FontWeight','bold')
title(['Date = ',time_utc,', Fraction (%) Solar Intensity Reduction on Sun Projected Disc due to an Obstruction of Radius = ',num2str(radius_obstruction),' km at SEL_1'], 'FontWeight','bold')
xlabel('Longitude [deg]'); ylabel('Latitude [deg]');
axis equal
xlim([-180 180]); ylim([-90 90]);

% Overlay 10° white grid
for lonV = lon_edges
    plot([lonV lonV], [-90 90], 'w', 'LineWidth', 0.5);
end
for latH = lat_edges
    plot([-180 180], [latH latH], 'w', 'LineWidth', 0.5);
end

% Annotate number of rays per cell (counts_total)

for r = 1:length(lat_centers)
    for c = 1:length(lon_centers)
        n = round(fraction_intensity_reduction_SRM(r,c),1);
        if n > 0 % Zeros are skipped
            ndust = counts_dust(r,c);
            text(lon_centers(c), lat_centers(r), sprintf('%s', num2str(n)), ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                'FontSize',7, 'FontWeight','bold', 'Color','w', ...
                'BackgroundColor',[0 0 0 0.4], 'Margin',1);
        end
    end
end

set(gca,'FontSize',14, 'FontName','Times New Roman')

% Save figure
saveas(fig3, replace(['Fraction_Intensity_Reduction_Projected_Disc_with_SRM_Obstruction_Rays',num2str(num_rays),'_Radius_',num2str(round(radius_obstruction)),'_Date_',time_utc,'.png'],':','_'));

