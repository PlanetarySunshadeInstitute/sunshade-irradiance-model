
function    [varargout] = ...
			........................................................................................................................................................ 
            analysis_General___0___X 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Start the program stopwatch, import the input model parameters and partitions, and begin a switch statement on the model settings.                      %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Start the program stopwatch.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


tic


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Import the input model and physical parameters.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


model_parameters = import_Model_Parameters___0___Sr;
[radii, ~, ~, ~] = parameters_Physical___0___2Sr2S;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	If the analysis looks just at irradiance, import the star and planet partitions. Otherwise import the shade partitions as well.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


if strcmp(model_parameters.analysis_type, 'irradiance')

	partitions   = generate_Partitions___A___Sr(model_parameters, {'star', 'planet'}, radii);

else 

	partitions   = generate_Partitions___A___Sr(model_parameters, {'star', 'planet', 'shade'}, radii);

end


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Begin a switch statement on the model settings.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Extract constellation metadata for the diagnostic chart.
%	shade partitions are only built when analysis_type is not 'irradiance'.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


if isfield(partitions, 'shade') && isfield(partitions.shade, 'heliogyros')
	num_shades = size(partitions.shade.heliogyros.initial.positions, 2);
else
	num_shades = [];
end


switch true


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Case: configuration 'matlab' (as the output is variable, the results will be displayed in the Matlab command window).                                   %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	case (strcmp(model_parameters.configuration, 'manual'))


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Identify a string to monitor program progress through the course of the analysis.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		counter_display = 'currently on time period: ';


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Perform the desired analysis.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		results         = analysis_Iterate_And_Switch___Sr2St___M (model_parameters.time, partitions, model_parameters.analysis_type, counter_display);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Output the results to the Matlab command window.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		varargout{1}    = results;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Compute and display the irradiance factor diagnostic chart.
%   A dedicated 12-step monthly analysis is run so the chart always spans the
%   full year regardless of how many time steps the main analysis used.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		fprintf('\nComputing diagnostic chart (12 monthly time steps)...\n');
		diag_year           = str2double(model_parameters.time.start(1:4));
		diag_time.start     = sprintf('%d-01-15T12:00:00', diag_year);
		diag_time.frequency = 30;
		diag_time.periods   = 12;

		diag_results        = analysis_Iterate_And_Switch___Sr2St___M ( ...
		                          diag_time, partitions, model_parameters.analysis_type, ...
		                          'diagnostic month: ');

		diag_lat_deg        = linspace(-90, 90, size(diag_results, 1));
		diag_central_col    = ceil(size(diag_results, 2) / 2);
		diag_results_2d     = reshape(diag_results(:, diag_central_col, :), ...
		                              size(diag_results, 1), size(diag_results, 3));
		diag_time_doy       = 15 + (0:11) * 30;   % mid-month: Jan 15, Feb 14, ...


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Apply solar declination coordinate correction to the diagnostic results.
%
%	The diagnostic uses the same heliocentric-to-geographic shift as the
%	preconfigured path, but is adapted for the manual-mode phi grid, which
%	is coarser (planet.phi.n cells over [-90, 90]).  The buffer depth is
%	computed dynamically so that it always covers the maximum declination
%	shift for the current grid resolution.
%
%	Note: with a coarse grid (e.g. 11 cells, step ~16.4 deg) the shift
%	toggles between -1, 0, and +1 cells across the year.  The diagnostic
%	chart will reflect this discrete stepping, which is the physically
%	correct behaviour at this resolution.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		n_phi_diag          = size(diag_results, 1);
		lat_step_diag       = 180 / n_phi_diag;                           % degrees per cell for this grid
		n_buf_diag          = ceil(23.44 / lat_step_diag) + 1;            % max shift + 1 cell margin

		diag_shifted_2d     = zeros(n_phi_diag, numel(diag_time_doy));

		for k = 1 : numel(diag_time_doy)
		    delta_k         = spencer_declination(diag_time_doy(k), 365);
		    shift_k         = round(delta_k / lat_step_diag);
		    buf_s_k         = repmat(diag_results_2d(1,   k), [n_buf_diag, 1]);
		    buf_n_k         = repmat(diag_results_2d(end, k), [n_buf_diag, 1]);
		    col_ext         = [buf_s_k; diag_results_2d(:, k); buf_n_k];
		    row_s_k         = n_buf_diag + 1 - shift_k;
		    diag_shifted_2d(:, k) = col_ext(row_s_k : row_s_k + n_phi_diag - 1);
		end

		excel_file          = location_Heliogyro_Kinematics_Data___E___Sr;
		plot_irradiance_diagnostic___M___0(diag_shifted_2d, diag_lat_deg, diag_time_doy, '', excel_file.name, num_shades);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Compute and display the area-weighted mean analysis type along the column,
%	reflecting the fact that Earth's rotation distributes shade across all longitudes.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%

		latitudes_deg   = linspace(-90, 90, size(results, 1));
		weights         = cosd(latitudes_deg);
		weights         = weights / sum(weights);
		central_col     = ceil(size(results, 2) / 2);
		central_strip   = results(:, central_col);
		mean_irradiance = sum(central_strip .* weights');
		% Diagnostic: show which constellation Excel input was used.
		excel_file       = location_Heliogyro_Kinematics_Data___E___Sr;
		fprintf('\nConstellation input: %s\n', excel_file.name);
		fprintf('\nArea-weighted mean %s: %.4f\n\n', model_parameters.analysis_type, mean_irradiance);


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Case: configuration 'preconfigured: low resolution' (results are automatically exported to the NC file, as they're too large for the command window).   %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	case (strcmp(model_parameters.configuration, 'preconfigured: low resolution'))


%%%%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%%%%
%%%%	Identify locations of the NC template and export files, and copy a template to be exported. 
%%%%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Identify the full locations of the NC template and export files.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		nc_files                        = location_NC_Templates_And_Exports___0___Sr;
		excel_file                      = location_Heliogyro_Kinematics_Data___E___Sr;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Copy and rename the template in the exports folder.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		copyfile(nc_files.locations.templates.irradiance_factor, nc_files.locations.exports.current_file);


%%%%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%%%%
%%%%	Analyze results for a regular year, and export them to the NC file. 
%%%%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Define the time start, frequency, and number of periods for a regular year.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		model_parameters.time.start     = model_parameters.time.regular_year.start;
		model_parameters.time.frequency = model_parameters.time.regular_year.frequency;
		model_parameters.time.periods   = model_parameters.time.regular_year.periods;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Identify a string to monitor program progress through the course of the analysis.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		counter_display                 = 'currently on regular year day: ';


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Analyze the results.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		results                         = analysis_Iterate_And_Switch___Sr2St___M (model_parameters.time, partitions, model_parameters.analysis_type, counter_display);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Apply solar declination coordinate correction (regular year).
%
%	The shade kernel computes irradiance factors in the heliocentric frame,
%	where phi = 0 always coincides with the sub-stellar point. To map those
%	values to the correct geographic latitude for each day, the result array
%	is shifted by round(delta / lat_step) rows, where delta is the solar
%	declination from Spencer (1971).
%
%	Buffer cells (n_buffer = 27 per hemisphere) are pre-filled with the value
%	of the outermost computed disc cell for each day.  This avoids a
%	discontinuity at the buffer boundary and is physically appropriate:
%	the buffer zone is consumed only near the solstices, when those latitudes
%	are in polar night and coszrs in CESM drives solin to zero regardless of
%	the sun_shade value.
%
%	Buffer depth (27 cells x 0.9424 deg/cell = 25.44 deg) exceeds the
%	maximum declination shift (23.44 deg / 0.9424 deg/cell = 24.87 cells,
%	rounded to 25) by 2 cells of margin.
%
%	Reference: Spencer, J.W. (1971). Search, 2(5), 172.
%	See also:  spencer_declination.m
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		n_buffer                        = 27;
		lat_step_deg                    = 0.942408376963357;              % degrees per latitude cell, matching the CESM f09 grid

		buffer_south                    = repmat(results(1,   :, :), [n_buffer, 1, 1]);
		buffer_north                    = repmat(results(end, :, :), [n_buffer, 1, 1]);
		results_extended                = [buffer_south; results; buffer_north];    % (192 + 2*n_buffer) x 1 x n_days

		n_days                          = size(results, 3);
		results_shifted                 = zeros(192, 1, n_days);

		for d = 1 : n_days
		    delta_deg                   = spencer_declination(d, 365);
		    shift                       = round(delta_deg / lat_step_deg);
		    row_s                       = n_buffer + 1 - shift;
		    assert(row_s >= 1 && row_s + 191 <= size(results_extended, 1), ...
		        'Declination shift exceeds buffer on regular-year day %d (shift = %d cells).', d, shift);
		    results_shifted(:, :, d)    = results_extended(row_s : row_s + 191, :, d);
		end

		clear results buffer_south buffer_north results_extended


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Extract 12 monthly samples from the declination-corrected results for the diagnostic chart.
%   This must happen before repmat/permute expand the array and before clear.
%   Indices correspond to approx. the 15th of each month in a 365-day year.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		diag_monthly_idx                = [15, 46, 74, 105, 135, 166, 196, 227, 257, 288, 318, 349];
		diag_lat_deg                    = linspace(-90, 90, size(results_shifted, 1));
		diag_results_2d                 = reshape(results_shifted(:, 1, diag_monthly_idx), ...
		                                          size(results_shifted, 1), numel(diag_monthly_idx));
		diag_time_doy                   = diag_monthly_idx;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Replicate the declination-corrected results across the second dimension (288 longitudinal coordinates), without replication (1) along the first and third dimensions.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		results_shifted                 = repmat(results_shifted, 1, 288, 1);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Permute the dimensions of the output matrix so that longitude, latitude, and time are represented along dimensions 1,2 and 3, respectively.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		results_shifted                 = permute(results_shifted, [2 1 3]);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Export results to the NC file.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		ncwrite(nc_files.locations.exports.current_file, 'df', results_shifted);
		clear results_shifted


%%%%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%%%%
%%%%	Analyze results for a leap year, and export them to the NC file.
%%%%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Define the time start, frequency, and number of periods for a leap year.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		model_parameters.time.start     = model_parameters.time.leap_year.start;
		model_parameters.time.frequency = model_parameters.time.leap_year.frequency;
		model_parameters.time.periods   = model_parameters.time.leap_year.periods;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Identify a string to monitor program progress through the course of the analysis.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		counter_display                 = 'currently on leap year day: ';


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Analyze the results.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		results                         = analysis_Iterate_And_Switch___Sr2St___M (model_parameters.time, partitions, model_parameters.analysis_type, counter_display);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Apply solar declination coordinate correction (leap year).
%
%	Same approach as the regular year.  days_in_year = 366 is passed to
%	spencer_declination so that the day-angle B is scaled to the correct
%	annual period for a leap year.
%
%	Note: the leap year run uses 364 periods (Jan 1 – Dec 30).  Day indices
%	1..364 are passed directly to spencer_declination; the formula remains
%	accurate to within its nominal ~0.035-degree tolerance across this range.
%
%	Reference: Spencer, J.W. (1971). Search, 2(5), 172.
%	See also:  spencer_declination.m
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		n_buffer                        = 27;
		lat_step_deg                    = 0.942408376963357;              % degrees per latitude cell, matching the CESM f09 grid

		buffer_south                    = repmat(results(1,   :, :), [n_buffer, 1, 1]);
		buffer_north                    = repmat(results(end, :, :), [n_buffer, 1, 1]);
		results_extended                = [buffer_south; results; buffer_north];    % (192 + 2*n_buffer) x 1 x n_days

		n_days                          = size(results, 3);
		results_shifted                 = zeros(192, 1, n_days);

		for d = 1 : n_days
		    delta_deg                   = spencer_declination(d, 366);
		    shift                       = round(delta_deg / lat_step_deg);
		    row_s                       = n_buffer + 1 - shift;
		    assert(row_s >= 1 && row_s + 191 <= size(results_extended, 1), ...
		        'Declination shift exceeds buffer on leap-year day %d (shift = %d cells).', d, shift);
		    results_shifted(:, :, d)    = results_extended(row_s : row_s + 191, :, d);
		end

		clear results buffer_south buffer_north results_extended


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Replicate the declination-corrected results across the second dimension (288 longitudinal coordinates), without replication (1) along the first and third dimensions.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		results_shifted                 = repmat(results_shifted, 1, 288, 1);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Permute the dimensions of the output matrix so that longitude, latitude, and time are represented along dimensions 1,2 and 3, respectively.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		results_shifted                 = permute(results_shifted, [2 1 3]);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Export results to the NC file.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		ncwrite(nc_files.locations.exports.current_file, 'dfl', results_shifted);
		clear results_shifted


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Rename the exported NC file to reflect the constellation type and export date.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		if     contains(excel_file.name, 'uniform',   'IgnoreCase', true), shading_code = 'U';
		elseif contains(excel_file.name, 'gaussian',  'IgnoreCase', true), shading_code = 'G';
		elseif contains(excel_file.name, 'polar',     'IgnoreCase', true), shading_code = 'P';
		elseif contains(excel_file.name, 'seasonal',  'IgnoreCase', true), shading_code = 'S';
		else,                                                                shading_code = 'X';
		end

		date_str                        = char(datetime('now', 'Format', 'yyyy-MM-dd'));
		[folder, ~, ~]                  = fileparts(nc_files.locations.exports.current_file);
		serial                          = 1;
		new_file                        = fullfile(folder, sprintf('L1-R%s-%s-%03d.nc', shading_code, date_str, serial));
		while isfile(new_file)
			serial                      = serial + 1;
			new_file                    = fullfile(folder, sprintf('L1-R%s-%s-%03d.nc', shading_code, date_str, serial));
		end
		movefile(nc_files.locations.exports.current_file, new_file);
		nc_files.locations.exports.current_file = new_file;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Save the irradiance factor diagnostic chart as a PNG alongside the NC file.
%   Uses the same filename stem as the NC file (e.g. L1-RU-2026-03-31-001.png).
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		[nc_folder, nc_stem, ~]         = fileparts(new_file);
		png_path                        = fullfile(char(nc_folder), [char(nc_stem), '.png']);
		plot_irradiance_diagnostic___M___0(diag_results_2d, diag_lat_deg, diag_time_doy, png_path, excel_file.name, num_shades);


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	End configuration cases and stop the program stopwatch.                                                                                                 %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	End output type cases.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


end


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	End the program stopwatch.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


toc


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










