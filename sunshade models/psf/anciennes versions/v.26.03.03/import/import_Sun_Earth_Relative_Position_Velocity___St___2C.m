
function    [position_SE_data, velocity_SE_data] = ...
			........................................................................................................................................................
            import_Sun_Earth_Relative_Position_Velocity___0___2C


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Import the sun-earth position and velocity vector associated with the input UTC date in the Excel and export them as a column vectors.                  %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Import the time defined in the function input_Model_Parameters___0___Sr.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


time_UTC                                                                         = input_Model_Parameters___0___Sr().time_UTC;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Import the Excel file location of the sun-earth position and velocity vs. time data and the corresponding data locations within the file.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


[excel_file_location, range_time_data, range_position_data, range_velocity_data] = location_and_ranges_Relative_Positions_Velocities_SE_Excel_File___0___3St;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Import the sun-earth position and velocity vs. time data.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


position_SE_data                                                                 = readmatrix(excel_file_location,'Range',range_position_data);
velocity_SE_data                                                                 = readmatrix(excel_file_location,'Range',range_velocity_data);
time_UTC_data                                                                    = readcell(excel_file_location,'Range',range_time_data);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Find the row corresponding to the input UTC date (represented by the first ten characters).     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


row_at_input_time_UTC                                                            = find(strncmp(time_UTC_data, time_UTC, 10));


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Collapse the sun-earth position and velocity data to the value at the input UTC date, and convert them a 3D column vectors.   
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


position_SE_data                                                                 = position_SE_data(row_at_input_time_UTC,:);
position_SE_data                                                                 = position_SE_data';

velocity_SE_data                                                                 = velocity_SE_data(row_at_input_time_UTC,:);
velocity_SE_data                                                                 = velocity_SE_data';


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%