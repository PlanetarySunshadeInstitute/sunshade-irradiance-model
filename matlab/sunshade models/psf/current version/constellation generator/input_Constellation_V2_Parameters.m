% input_Constellation_V2_Parameters.m
%
% USER SETTINGS for the V2 time-varying constellation generator.
% Edit this file, then call generate_constellation_v2.m.
%
% ARCHITECTURE
%   A constellation is made up of one or more 'blobs'. Each blob is a
%   group of heliogyros sharing an elliptical footprint in the Y-Z plane,
%   stacked across multiple X-depth planes. The blob moves as a rigid body
%   through SEL1 space over the course of the year, defined by a motion
%   function you supply.
%
% BLOB PARAMETERS
%   center_position  [X, Y, Z] in km (L1-centred synodic frame).
%                    Default X = 813400 km (sunward of SEL1).
%                    Y and Z set the nominal ellipse center in the Y-Z plane.
%
%   ellipse_radii    [Ry, Rz] in km — semi-axes of the blob footprint in
%                    the Y-Z plane.
%
%   n_craft          Total spacecraft in this blob.
%
%   n_planes         Number of X-depth planes.
%
%   plane_spacing_km Spacing between adjacent planes (km). Spacing must
%                    exceed min_buffer_km to guarantee cross-plane separation.
%
%   min_buffer_km    Minimum center-to-center spacing between craft within
%                    the same plane (km). A warning is raised if this cannot
%                    be satisfied after the placement attempts.
%
%   motion_function  Function handle @(t) [dX, dY, dZ] — displacement of
%                    the blob center from center_position at day t (1-indexed,
%                    where t=1 is January 1). Return value is a 1x3 or 3x1
%                    vector in km. Use @(t) [0,0,0] for a stationary blob.
%
% MOTION FUNCTION EXAMPLES
%   Stationary:
%       @(t) [0, 0, 0]
%
%   Annual Z oscillation (e.g. seasonal north/south sweep of ±A km):
%       @(t) [0, 0, A * sin(2*pi*(t-1)/365)]
%
%   Phase-shifted oscillation (e.g. peaks at summer solstice, day ~172):
%       @(t) [0, 0, A * sin(2*pi*(t-172)/365)]
%
% COORDINATE CONVENTION
%   X positive = sunward.  Y = tangential.  Z = ecliptic north (out-of-plane).
%   All positions in km, L1-centred synodic frame.
%
% Authors: Planetary Sunshade Foundation

function params = input_Constellation_V2_Parameters()

% =========================================================================
% GLOBAL SETTINGS
% =========================================================================

params.year_type      = 'normal';   % 'normal' (365 days) or 'leap' (366 days)
params.label          = 'B2';       % Short label used in the output filename
                                    % e.g. 'B2' -> constellation_B2_11000crafts_..._001.mat

% =========================================================================
% BLOB DEFINITIONS
% =========================================================================
% Add or remove blobs by copying/deleting the block below.
% Each blob is params.blobs(n) — a struct with the fields described above.

% --- Blob 1: North polar blob ---
params.blobs(1).center_position  = [813400, 0, 840];  % km, X=sunward of SEL1
params.blobs(1).ellipse_radii    = [3000, 1800];        % [Ry, Rz] semi-axes, km
params.blobs(1).n_craft          = 7365;
params.blobs(1).n_planes         = 6;
params.blobs(1).plane_spacing_km = 10000;
params.blobs(1).min_buffer_km    = 45;
params.blobs(1).motion_function  = @(t) [0, 0, 3850*sin(2*pi*t/365 - 0.423*pi)];    

% --- Blob 2: South polar blob ---
params.blobs(2).center_position  = [813400, 0, 17]; % km
params.blobs(2).ellipse_radii    = [2000, 1200];        % [Ry, Rz] semi-axes, km
params.blobs(2).n_craft          = 3191;
params.blobs(2).n_planes         = 6;
params.blobs(2).plane_spacing_km = 10000;
params.blobs(2).min_buffer_km    = 45;
params.blobs(2).motion_function  = @(t) [0, 0, 11100*sin(2*pi*t/365 - 0.4*pi)]

% --- Blob 3: equatorial blob ---
params.blobs(3).center_position  = [813400, 0, -100]; % km
params.blobs(3).ellipse_radii    = [3000, 1800];        % [Ry, Rz] semi-axes, km
params.blobs(3).n_craft          = 0;
params.blobs(3).n_planes         = 6;
params.blobs(3).plane_spacing_km = 10000;
params.blobs(3).min_buffer_km    = 45;
params.blobs(3).motion_function  = @(t) [0, 0, 0];     

end
