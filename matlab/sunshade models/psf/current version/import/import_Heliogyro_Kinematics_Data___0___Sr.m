
function    [heliogyros] = ...
			........................................................................................................................................................
            import_Heliogyro_Kinematics_Data___0___Sr


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Import the heliogyro kinematics data.                                                                                       	                        %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Import the file location of the Excel file containing the heliogyro kinematics data.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


excel_file                                                                  = location_Heliogyro_Kinematics_Data___E___Sr;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Import the heliogyro initial positions and normal vectors.
%	Branch on shade file format: 'mat_v2' loads from a V2 .mat file; 'xlsx_v1' uses the original Excel path.
%	For mat_v2, day 1 is loaded here as a default. analysis_Iterate_And_Switch refreshes per-day during the run.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


if isfield(excel_file, 'format') && strcmp(excel_file.format, 'mat_v2')


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	V2 .mat constellation: load positions and normals for day 1 as initial values.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


	mat_data                                                                = load(excel_file.location);
	heliogyros.initial.positions                                            = mat_data.positions(:, :, 1);
	heliogyros.initial.normal_vectors                                       = mat_data.normals(:, :, 1);
	heliogyros.initial.time                                                 = {''};


else


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	V1 .xlsx constellation: read position vectors, normal vectors, and initial time from Excel ranges.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


	heliogyros.initial.positions                                            = readmatrix(excel_file.location,'Range', excel_file.ranges.position_vectors)';
	heliogyros.initial.normal_vectors                                       = readmatrix(excel_file.location,'Range', excel_file.ranges.normal_vectors)';
	heliogyros.initial.time                                                 = readcell(excel_file.location,'Range', excel_file.ranges.time_initial);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Remove all entries imported as NaN from the heliogyro position and normal vector sets.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


	heliogyros.initial.positions                                            = heliogyros.initial.positions(:, ~any(isnan(heliogyros.initial.positions)));
	heliogyros.initial.normal_vectors                                       = heliogyros.initial.normal_vectors(:, ~any(isnan(heliogyros.initial.normal_vectors)));


end


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	As the heliogryos are flat surfaces, choose their normals to be represented oriented towards the Star (for computational regularity).     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


heliogyros.initial.normal_vectors(:, heliogyros.initial.normal_vectors(1,:)<0) = -heliogyros.initial.normal_vectors(:, heliogyros.initial.normal_vectors(1,:)<1);


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%