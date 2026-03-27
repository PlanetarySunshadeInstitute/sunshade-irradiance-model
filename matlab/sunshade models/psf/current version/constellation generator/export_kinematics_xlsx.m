% export_kinematics_xlsx.m
%
% Exports constellation spacecraft positions to the heliogyro kinematics
% Excel file format expected by the PSF irradiance model.
%
% Writes to: excel/psf model/constellation_[pattern]_[N]craft_[date].xlsx
%
% The output format matches the existing heliogyro kinematics data.xlsx
% schema exactly:
%   PX, PY, PZ : position in km, L1-centred synodic frame
%   NX, NY, NZ : surface normal unit vector (all [1,0,0] for now)
%
% USAGE
%   filepath = export_kinematics_xlsx(positions, params, excel_folder)
%
% INPUT
%   positions    : struct from place_spacecraft.m
%   params       : struct from define_constellation_params.m
%   excel_folder : path to the excel/psf model/ directory
%
% OUTPUT
%   filepath     : full path to the written file
%
% Authors: Planetary Sunshade Foundation

function filepath = export_kinematics_xlsx(positions, params, excel_folder)

% -------------------------------------------------------------------------
% BUILD FILENAME
% -------------------------------------------------------------------------
if isfield(params, 'output_filename') && ~isempty(params.output_filename)
    filename = params.output_filename;
else
    date_str  = datestr(now, 'yyyy-mm-dd');
    filename  = sprintf('constellation_%s_%dcrafts_%s.xlsx', ...
                         params.pattern, params.N, date_str);
end
filepath  = fullfile(excel_folder, filename);

fprintf('[export] Writing %d spacecraft to:\n  %s\n', params.N, filepath);

% -------------------------------------------------------------------------
% BUILD TABLE
% -------------------------------------------------------------------------
% Header row matches existing heliogyro kinematics data.xlsx schema exactly.
% Rows 1-10 are reserved for metadata (read by the irradiance model as
% blank / ignored). Data starts at row 11, consistent with the existing
% file location_Heliogyro_Kinematics_Data___E___Sr.m range 'A11:C10000'.

header = {'PX', 'PY', 'PZ', 'NX', 'NY', 'NZ'};

data = [positions.PX, positions.PY, positions.PZ, ...
        positions.NX, positions.NY, positions.NZ];

T = array2table(data, 'VariableNames', header);

% -------------------------------------------------------------------------
% WRITE METADATA ROWS (rows 2-10)
% -------------------------------------------------------------------------
% These rows sit above the data and are ignored by the irradiance model's
% range reader. They document the constellation for human readers of the
% Excel file.

if strcmp(params.pattern, 'uniform')
    footprint_line = sprintf('Footprint:  %.0f km radius  (%.0f%% of solar disk)', ...
                             params.constellation_radius_km, ...
                             params.footprint_pct_of_disk);
else
    footprint_line = sprintf('Footprint:  polar ellipses Y %.0f km, Z %.0f km', ...
                             params.polar_radius_y_km, params.polar_radius_z_km);
end

meta = { ...
    'PSF Constellation Generator — heliogyro kinematics data', '', '', '', '', '';
    sprintf('Generated:  %s', datestr(now, 'yyyy-mm-dd HH:MM:SS')),  '', '', '', '', '';
    sprintf('Pattern:    %s', params.pattern),                        '', '', '', '', '';
    sprintf('N craft:    %d', params.N),                              '', '', '', '', '';
    sprintf('N planes:   %d  (spacing: %.0f km)', ...
             params.n_planes, params.plane_spacing_km),               '', '', '', '', '';
    footprint_line,                                                   '', '', '', '', '';
    sprintf('Profile:    %s', params.footprint_profile),              '', '', '', '', '';
    sprintf('Sail radius: %.0f km  (min buffer: %.0f km)', ...
             params.sail_radius_km, params.min_buffer_km),            '', '', '', '', '';
    sprintf('Coord:      L1-centred synodic km  (X+ sunward)'),       '', '', '', '', '';
};

% Write metadata to rows 1-9, header to row 10, data from row 11
% We use writecell for the metadata block and writetable for the data,
% writing to the same sheet sequentially.

if exist(filepath, 'file')
    delete(filepath);
end

% Write metadata block (rows 1-9)
writecell(meta, filepath, ...
          'Sheet', 1, ...
          'Range', 'A1');

% Write column headers (row 10)
writecell(header, filepath, ...
          'Sheet', 1, ...
          'Range', 'A10');

% Write data (rows 11 onward)
writematrix(data, filepath, ...
            'Sheet', 1, ...
            'Range', 'A11');

% -------------------------------------------------------------------------
% VERIFY
% -------------------------------------------------------------------------
if exist(filepath, 'file')
    info     = dir(filepath);
    size_kb  = info.bytes / 1024;
    fprintf('[export] Done. File size: %.1f KB\n', size_kb);
    fprintf('[export] Rows written: %d (data) + 10 (header/metadata)\n', params.N);
else
    error('[export] File not found after write — check excel_folder path:\n  %s', ...
          excel_folder);
end

% -------------------------------------------------------------------------
% POLAR METADATA (appended if relevant)
% -------------------------------------------------------------------------
if strcmp(params.pattern, 'polar')
    polar_meta = { ...
        '', '', '', '', '', '';
        sprintf('Polar targeting:'), '', '', '', '', '';
        sprintf('  Polar center Z:    %.1f km (%s)', ...
                 params.polar_center_z_km, params.hemisphere), '', '', '', '', '';
        sprintf('  Polar fraction:    %.0f%%', ...
                 params.polar_fraction * 100),                   '', '', '', '', '';
        sprintf('  Z center offset:   %.1f km', ...
                 params.polar_center_z_km),                      '', '', '', '', '';
        sprintf('  Ellipse axes:      Y %.1f km, Z %.1f km', ...
                 params.polar_radius_y_km, ...
                 params.polar_radius_z_km),                      '', '', '', '', '';
    };

    % Append after data block
    next_row  = params.N + 12;
    range_str = sprintf('A%d', next_row);
    writecell(polar_meta, filepath, 'Sheet', 1, 'Range', range_str);
end

fprintf('\n');

end