function paths = config_paths()
%CONFIG_PATHS  Central path configuration for the sunshade-irradiance project.
%
%   paths = config_paths() returns a struct of absolute folder paths,
%   all derived from a single project root so the codebase runs on any
%   machine without editing individual scripts.
%
% HOW IT WORKS
%   1. Checks the environment variable SUNSHADE_PROJECT_ROOT.
%      Set this once per machine (e.g. in startup.m or your shell profile)
%      if the repo is not at a standard location.
%
%   2. If the variable is absent, falls back to the folder containing
%      this file (config_paths.m lives at the repo root, so this gives
%      the correct answer automatically after a fresh clone).
%
% FIRST-TIME SETUP
%   Option A — no setup needed:
%     Clone the repo, add it to the MATLAB path.  config_paths() will
%     auto-detect the root from its own location.
%
%   Option B — explicit root (e.g. if you have multiple clones):
%     In your MATLAB startup.m or shell profile, set:
%       setenv('SUNSHADE_PROJECT_ROOT', '/path/to/your/clone/MatLab')
%
% PATHS RETURNED
%   paths.root             — repo root (all others are derived from this)
%   paths.excel_folder     — excel/psf model/
%   paths.nc_folder        — numerical control/
%   paths.nc_exports       — numerical control/exports/
%   paths.constellation    — matlab/.../constellation generator/
%   paths.file_locations   — matlab/.../file locations/
%
% Authors: Planetary Sunshade Foundation


    % ---- 1. Resolve project root ----------------------------------------
    root = getenv('SUNSHADE_PROJECT_ROOT');

    if isempty(root)
        % Fallback: derive from this file's own location.
        % config_paths.m lives at the repo root, so fileparts gives the root.
        root = fileparts(mfilename('fullpath'));
    end

    % Normalise: strip any trailing slash so fullfile works uniformly.
    root = strtrim(root);
    if ~isempty(root) && (root(end) == '/' || root(end) == '\')
        root = root(1:end-1);
    end

    if ~isfolder(root)
        error(['config_paths: project root does not exist: %s\n' ...
               'Set the SUNSHADE_PROJECT_ROOT environment variable to the ' ...
               'absolute path of your MatLab repo clone.'], root);
    end

    % ---- 2. Build derived paths -----------------------------------------
    paths.root           = root;
    paths.excel_folder   = fullfile(root, 'excel',    'psf model');
    paths.nc_folder      = fullfile(root, 'numerical control');
    paths.nc_exports     = fullfile(root, 'numerical control', 'exports');
    paths.constellation  = fullfile(root, 'matlab', 'sunshade models', ...
                                    'psf', 'current version', 'constellation generator');
    paths.file_locations = fullfile(root, 'matlab', 'sunshade models', ...
                                    'psf', 'current version', 'file locations');

end
