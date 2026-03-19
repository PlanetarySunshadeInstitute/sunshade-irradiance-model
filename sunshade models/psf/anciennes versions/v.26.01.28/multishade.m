% Multi-Object Penumbral Shading Simulation on Earth  (3-axis offsets + progress)
% Computes and visualizes cumulative penumbral shading (%) on Earth's
% surface due to many annular objects, reporting percent-complete status,
% and displays average sunlight reduction over the Earth disk.
%

clear; clc;

%% ------------------------ 1.  INPUTS ------------------------------------
inner_radius = 10;          % annulus inner radius (km)
outer_radius = 20;          % annulus outer radius (km)
num_objs     = 8;         % number of shading objects

% ── Distance distribution (10^6 km) ──────────────────────────────────────
% Examples – use any you like; a few common patterns are shown below:
%   Uniform linear:     track_func = @(N) linspace(track_min,track_max,N);
%   Random uniform:     track_func = @(N) track_min + (track_max-track_min).*rand(1,N);
%   Cosine U-shape:     track_func = @(N) track_min + (track_max-track_min)*((1-nthroot(cos(linspace(0,pi,N)),1))/2);

D_min       = 2.31;
D_max       = 2.41;
dist_func = @(N) D_min + (D_max-D_min).*rand(1,N);
D_mill_km = dist_func(num_objs);

% ── Cross-track / "normal" offset (Δy, km) ───────────────────────────────
off_min     = -13000;
off_max     = 13000;
offset_func = @(N) off_min + (off_max-off_min)*((1 - (nthroot(cos((linspace(0,pi,N)-0)),3))) / 2);
offset_km = offset_func(num_objs);

% ── Down-track / "along-orbit" offset (Δx, km) ──────────────────────
track_min   = -6000;
track_max   = 6000;
track_func = @(N) track_min + (track_max-track_min).*rand(1,N);
track_km = track_func(num_objs);

%% ------------------------ 2.  CONSTANTS & GRID --------------------------
R_sun_km  = 695700;
AU_km     = 149597870;
R_sun_ang = asin(R_sun_km / AU_km);        % Sun's angular radius (rad)
EarthRadius  = 6371;

% Ground grid large enough to capture full penumbra of furthest object
n      = 1000;                                           % grid resolution
r_max_penumbra = (min(D_mill_km) * 1e6) * (R_sun_ang + outer_radius / (min(D_mill_km) * 1e6));
r_max  = max(r_max_penumbra, EarthRadius);
x      = linspace(-r_max, r_max, n);
y      = linspace(-r_max, r_max, n);
[X,Y]  = meshgrid(x,y);

shading_total = zeros(n,n);
pctStep       = max(floor(num_objs/100),1);             % 1 % reporting step

%% ------------------------ 3.  OVERLAP AREA HELPER -----------------------
function A = overlap_area(R, r, d)
    if d >= (R + r)
        A = 0;
    elseif d <= abs(R - r)
        A = pi * min(R, r)^2;
    else
        alpha = acos((d.^2 + R.^2 - r.^2) ./ (2 .* d .* R));
        beta  = acos((d.^2 + r.^2 - R.^2) ./ (2 .* d .* r));
        A = R.^2 .* alpha + r.^2 .* beta ...
            - 0.5 .* sqrt(( -d + R + r ) .* ( d + R - r ) ...
            .* ( d - R + r ) .* ( d + R + r ));
    end
end


%% ------------------------ 4.  PROGRESS DISPLAY --------------------------
dq = parallel.pool.DataQueue;
afterEach(dq,@printPct);
function printPct(p)
    persistent last; if isempty(last), last = 0; end
    if p > last, last = p; fprintf('Progress: %d%% complete\n',p); end
end

%% ------------------------ 5.  MAIN LOOP (PARFOR) ------------------------
parfor k = 1:num_objs
    % Progress message
    if mod(k,pctStep)==0 || k==num_objs, send(dq,round(100*k/num_objs)); end

    % Current object's parameters
    D_km   = D_mill_km(k)*1e6;
    dy_km  = offset_km(k);
    dx_km  = track_km(k);

    Rout   = outer_radius / D_km;          % outer annulus angular radius
    Rin    = inner_radius / D_km;          % inner annulus angular radius

    % Angular separation for each ground point (3-axis)
    Sep_ang = sqrt( (X+dx_km).^2 + (Y+dy_km).^2 ) / D_km;

    % Overlap fractions
    over_out = arrayfun(@(d) overlap_area(Rout,R_sun_ang,d), Sep_ang);
    over_in  = arrayfun(@(d) overlap_area(Rin ,R_sun_ang,d), Sep_ang);
    shade_i  = (over_out - over_in) ./ (pi*R_sun_ang^2) * 100;

    % Variable opacity (inner=1, outer= inner/outer)
    t       = (Sep_ang - Rin) ./ (Rout - Rin);
    t       = min(max(t,0),1);
    opacity = 1 - (1 - (inner_radius/outer_radius)) .* t;
    shade_i = shade_i .* opacity;

    shade_i(Sep_ang > (Rout + R_sun_ang)) = 0;  % outside penumbra
    shading_total = shading_total + shade_i;    % accumulate
end

%% ------------------------ 6.  METRICS & VISUALS -------------------------

mask_earth   = (X.^2 + Y.^2) <= EarthRadius^2;
avg_reduction= mean(shading_total(mask_earth),'all')

