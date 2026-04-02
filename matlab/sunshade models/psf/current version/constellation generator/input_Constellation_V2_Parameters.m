% input_Constellation_V2_Parameters.m
%
% USER SETTINGS for the V2 time-varying constellation generator.
% Edit this file, then call generate_constellation_v2.m.
%
% ARCHITECTURE
%   A constellation is made up of one or more 'clusters'. Each cluster is
%   a group of heliogyros sharing an elliptical footprint in the Y-Z plane.
%   The cluster moves as a rigid body through SEL1 space over the year,
%   defined by a motion function you supply.
%
%   Clusters are named (e.g. 'Blue Cluster') and their diagnostic chart
%   colors correspond to MATLAB's lines() colormap in definition order:
%     Cluster 1 → blue,  Cluster 2 → red,  Cluster 3 → yellow, etc.
%   Set active = false to exclude a cluster from the run without removing
%   its definition — equivalent to the old n_craft = 0 convention.
%
% PLACEMENT METHOD
%   'lattice'  — Hexagonal close-packing with FCC 3D offset across planes.
%                Craft count per cluster is determined by geometry.
%                No buffer violation checks are needed — analytically guaranteed.
%
%   'random'   — Random sampling within the ellipse with minimum buffer
%                enforcement. Craft count per cluster is set by n_craft.
%
% CLUSTER PARAMETERS (all methods)
%   name             Display name shown in diagnostic charts and console.
%   active           true = include in run.  false = skip (zero craft).
%   center_position  [X, Y, Z] in km (L1-centred synodic frame).
%                    X = 838400 km is the optimal shading distance.
%   ellipse_radii    [Ry, Rz] semi-axes of the cluster footprint in Y-Z (km).
%   motion_function  @(t) [dX, dY, dZ] — displacement from center_position
%                    at day t (1-indexed, t=1 = Jan 1). Returns 1x3 or 3x1.
%
% CLUSTER PARAMETERS (random method only)
%   n_craft          Total spacecraft to place in this cluster.
%
% MOTION FUNCTION EXAMPLES
%   Stationary:
%       @(t) [0, 0, 0]
%   Annual Z oscillation (±A km):
%       @(t) [0, 0, A * sin(2*pi*(t-1)/365)]
%   Phase-shifted (peaks at summer solstice, day ~172):
%       @(t) [0, 0, A * sin(2*pi*(t-172)/365)]
%
% COORDINATE CONVENTION
%   X positive = sunward.  Y = tangential.  Z = ecliptic north.
%   All positions in km, L1-centred synodic frame.
%
% Authors: Planetary Sunshade Foundation

function params = input_Constellation_V2_Parameters()

% =========================================================================
% GLOBAL SETTINGS
% =========================================================================

params.year_type = 'normal';   % 'normal' (365 days) or 'leap' (366 days)
params.label     = 'B2';       % Short label used in the output filename

% =========================================================================
% PLACEMENT METHOD AND PARAMETERS
% =========================================================================

params.placement.method = 'lattice';   % 'lattice' (default) or 'random'

% --- Lattice parameters (used when method = 'lattice') ---
%
%   n_planes          Number of X-depth planes. Sweet spot is 3 — gives
%                     one FCC layer per plane (A, B, C) with zero line-of-
%                     sight overlap. A 4th plane repeats layer A.
%
%   plane_spacing_km  Distance between adjacent planes (km). Any value
%                     from 5 to ~500 km preserves zero shadow overlap.
%
%   min_buffer_km     Edge-to-edge clearance between shades within a plane.
%                     Hex spacing = 2*sail_radius + min_buffer.
%
%   sail_radius_km    Physical radius of each heliogyro (km).
%                     Must match input_Model_Parameters.
%
%   los_margin_km     Extra cross-layer clearance margin beyond 2*sail_radius.
%                     Default 0 is geometrically sufficient.

params.placement.lattice.n_planes         = 3;
params.placement.lattice.plane_spacing_km = 500;
params.placement.lattice.min_buffer_km    = 45;
params.placement.lattice.sail_radius_km   = 20;
params.placement.lattice.los_margin_km    = 0;

% --- Random parameters (used when method = 'random') ---
%   min_buffer_km is center-to-center here (not edge-to-edge).

params.placement.random.n_planes         = 6;
params.placement.random.plane_spacing_km = 10000;
params.placement.random.min_buffer_km    = 45;

% =========================================================================
% CLUSTER DEFINITIONS
% =========================================================================
% Add or remove clusters by copying/deleting the blocks below.
% Set active = false to exclude a cluster without deleting its definition.
% Colors in the diagnostic correspond to MATLAB lines() order: blue, red, yellow, ...

% --- Cluster 1: Blue Cluster ---
params.clusters(1).name             = 'Blue Cluster';
params.clusters(1).active           = true;
params.clusters(1).center_position  = [813400, 0,  840];
params.clusters(1).ellipse_radii    = [2200, 1500];
params.clusters(1).motion_function  = @(t) [0, 0, 3850*sin(2*pi*t/365 - 0.423*pi)];
params.clusters(1).n_craft          = 7365;   % random method only

% --- Cluster 2: Red Cluster ---
params.clusters(2).name             = 'Red Cluster';
params.clusters(2).active           = true;
params.clusters(2).center_position  = [813400, 0,   17];
params.clusters(2).ellipse_radii    = [1500, 1000];
params.clusters(2).motion_function  = @(t) [0, 0, 11100*sin(2*pi*t/365 - 0.4*pi)];
params.clusters(2).n_craft          = 3191;   % random method only

% --- Cluster 3: Yellow Cluster ---
params.clusters(3).name             = 'Yellow Cluster';
params.clusters(3).active           = false;
params.clusters(3).center_position  = [813400, 0, -100];
params.clusters(3).ellipse_radii    = [3000, 1800];
params.clusters(3).motion_function  = @(t) [0, 0, 0];
params.clusters(3).n_craft          = 0;      % random method only

end
