
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
%	Import the heliogyro initial positions and normal vectors as separate matrices, and the initial time as a string.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


heliogyros.initial.positions                                                = readmatrix(excel_file.location,'Range', excel_file.ranges.position_vectors)';
heliogyros.initial.normal_vectors                                           = readmatrix(excel_file.location,'Range', excel_file.ranges.normal_vectors)';
heliogyros.initial.time                                                     = readcell(excel_file.location,'Range', excel_file.ranges.time_initial);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Remove all entries imported as NaN from the heliogyro position and normal vector sets.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


heliogyros.initial.positions                                                = heliogyros.initial.positions(:, ~any(isnan(heliogyros.initial.positions)));
heliogyros.initial.normal_vectors                                           = heliogyros.initial.normal_vectors(:, ~any(isnan(heliogyros.initial.normal_vectors)));


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	As the heliogryos are flat surfaces, choose their normals to be represented oriented towards the Star (for computational regularity).     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


heliogyros.initial.normal_vectors(:, heliogyros.initial.normal_vectors(1,:)<0) = -heliogyros.initial.normal_vectors(:, heliogyros.initial.normal_vectors(1,:)<1);


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%