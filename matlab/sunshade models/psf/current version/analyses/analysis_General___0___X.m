
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


		excel_file          = location_Heliogyro_Kinematics_Data___E___Sr;
		plot_irradiance_diagnostic___M___0(diag_results_2d, diag_lat_deg, diag_time_doy, '', excel_file.name, num_shades);


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


		fprintf('DEBUG: Template path = %s\n', nc_files.locations.templates.irradiance_factor);
		fprintf('DEBUG: Export path = %s\n', nc_files.locations.exports.current_file);
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
%	Extract 12 monthly samples for the diagnostic chart.
%   Indices correspond to approx. the 15th of each month in a 365-day year.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		diag_monthly_idx                = [15, 46, 74, 105, 135, 166, 196, 227, 257, 288, 318, 349];
		diag_lat_deg                    = linspace(-90, 90, size(results, 1));
		diag_results_2d                 = reshape(results(:, 1, diag_monthly_idx), ...
		                                          size(results, 1), numel(diag_monthly_idx));
		diag_time_doy                   = diag_monthly_idx;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Replicate results across the second dimension (288 longitudinal coordinates), without replication (1) along the first and third dimensions.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		results                         = repmat(results, 1, 288, 1);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Permute the dimensions of the output matrix so that longitude, latitude, and time are represented along dimensions 1,2 and 3, respectively.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		results                         = permute(results, [2 1 3]);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Export results to the NC file.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		ncwrite(nc_files.locations.exports.current_file, 'df', results);
		clear results


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
%	Replicate results across the second dimension (288 longitudinal coordinates), without replication (1) along the first and third dimensions.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		results                         = repmat(results, 1, 288, 1);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Permute the dimensions of the output matrix so that longitude, latitude, and time are represented along dimensions 1,2 and 3, respectively.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		results                         = permute(results, [2 1 3]);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Export results to the NC file.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		ncwrite(nc_files.locations.exports.current_file, 'dfl', results);
		clear results


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
%   The vertical axis shows degrees from the sub-stellar point (sun-frame phi),
%   not geographic latitude — CESM handles Earth rotation separately.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		[nc_folder, nc_stem, ~]         = fileparts(new_file);
		png_path                        = fullfile(char(nc_folder), [char(nc_stem), '.png']);
		plot_irradiance_diagnostic___M___0(diag_results_2d, diag_lat_deg, diag_time_doy, png_path, excel_file.name, num_shades);


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Case: configuration 'preconfigured: full resolution'                                                                                                    %&%%
%%%%	Computes irradiance factors across all 144 disc longitude columns (the full sun-facing hemisphere), producing a physically correct 192 x 288 output.    %&%%
%%%%	The sub-stellar point is placed at CESM longitude 0°. The 144 nightside columns (lon 90°–270°) are filled with 1.0 (no shading).                        %&%%
%%%%	CESM rotates the Earth underneath this fixed sun-frame pattern when applying the file in a climate run.                                                  %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	case (strcmp(model_parameters.configuration, 'preconfigured: full resolution'))


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
%   Output is 192 (lat) x 144 (disc theta) x 365 (days).
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		results                         = analysis_Iterate_And_Switch___Sr2St___M (model_parameters.time, partitions, model_parameters.analysis_type, counter_display);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Extract 12 monthly samples for the diagnostic chart using the central disc column (theta ≈ 0°, the sub-stellar meridian).
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		diag_monthly_idx                = [15, 46, 74, 105, 135, 166, 196, 227, 257, 288, 318, 349];
		diag_lat_deg                    = linspace(-90, 90, size(results, 1));
		diag_central_col                = ceil(size(results, 2) / 2);
		diag_results_2d                 = reshape(results(:, diag_central_col, diag_monthly_idx), ...
		                                          size(results, 1), numel(diag_monthly_idx));
		diag_time_doy                   = diag_monthly_idx;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Map the 192 x 144 disc result into the 192 x 288 CESM longitude grid.
%
%   Convention: sub-stellar point at CESM longitude 0° (computation at noon UTC).
%
%   Disc bins 73–144 (theta = +0.625° to +89.375°) → CESM lon bins   1–72  (0°    to  88.75°, dayside east)
%   Disc bins  1–72  (theta = -89.375° to -0.625°) → CESM lon bins 217–288 (270°  to 358.75°, dayside west)
%   CESM lon bins 73–216 (90° to 268.75°)          → 1.0                            (nightside, no shading)
%
%   Corner cells where sin²(theta) + sin²(phi) > 1 are outside the Earth's disc edge and the
%   irradiance computation already returns 1.0 for them — no additional masking is needed.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		full_results                                = ones(size(results,1), 288, size(results,3));
		full_results(:,   1:72,  :)                 = results(:, 73:144, :);   % eastern dayside
		full_results(:, 217:288, :)                 = results(:,  1:72,  :);   % western dayside
		results                                     = full_results;
		clear full_results


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Permute the dimensions of the output matrix so that longitude, latitude, and time are represented along dimensions 1, 2, and 3, respectively.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		results                         = permute(results, [2 1 3]);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Export results to the NC file.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		ncwrite(nc_files.locations.exports.current_file, 'df', results);
		clear results


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
%   Output is 192 (lat) x 144 (disc theta) x 366 (days).
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		results                         = analysis_Iterate_And_Switch___Sr2St___M (model_parameters.time, partitions, model_parameters.analysis_type, counter_display);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Map the 192 x 144 disc result into the 192 x 288 CESM longitude grid (same convention as regular year above).
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		full_results                                = ones(size(results,1), 288, size(results,3));
		full_results(:,   1:72,  :)                 = results(:, 73:144, :);   % eastern dayside
		full_results(:, 217:288, :)                 = results(:,  1:72,  :);   % western dayside
		results                                     = full_results;
		clear full_results


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Permute the dimensions of the output matrix so that longitude, latitude, and time are represented along dimensions 1, 2, and 3, respectively.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		results                         = permute(results, [2 1 3]);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Export results to the NC file.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


		ncwrite(nc_files.locations.exports.current_file, 'dfl', results);
		clear results


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
%   Uses the central disc column (sub-stellar meridian) for the sub-stellar degree profile.
%   The vertical axis shows degrees from the sub-stellar point (sun-frame phi),
%   not geographic latitude — CESM handles Earth rotation separately.
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










