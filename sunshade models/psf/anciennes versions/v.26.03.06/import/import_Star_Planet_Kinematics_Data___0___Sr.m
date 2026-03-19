
function    [kinematics_SP_data] = ...
			........................................................................................................................................................
            import_Star_Planet_Kinematics_Data___0___Sr


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Import the star planet kinematics data at the input UTC datem perform calculations, and export all results in a structure.                              %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Import the time defined in the function input_Model_Parameters___0___Sr.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


time_UTC                        = input_Model_Parameters___0___Sr().time_UTC;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Import the Excel file location of the star planet kinematics data vs. time data and the corresponding data locations within the file.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


excel_file                      = location_Star_Planet_Kinematics_Data___E___Sr;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Import the star planet kinematics vs. time data.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


kinematics_SP_data.position     = readmatrix(excel_file.location,'Range', excel_file.ranges.position_data);
kinematics_SP_data.velocity     = readmatrix(excel_file.location,'Range', excel_file.ranges.velocity_data);
kinematics_SP_data.time_UTC     = readcell(excel_file.location,'Range', excel_file.ranges.time_data);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Find the row corresponding to the input UTC date (represented by the first ten characters).     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


row_at_input_time_UTC           = find(strncmp(kinematics_SP_data.time_UTC, time_UTC, 10));


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Collapse the star earth kinematics data to the value at the input UTC date, and convert them a 3D column vectors.   
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


kinematics_SP_data.position     = kinematics_SP_data.position(row_at_input_time_UTC,:);
kinematics_SP_data.position     = kinematics_SP_data.position';

kinematics_SP_data.velocity     = kinematics_SP_data.velocity(row_at_input_time_UTC,:);
kinematics_SP_data.velocity     = kinematics_SP_data.velocity';


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the star planet distance and the distance to L1 (the latter is added to the structure within the called function).      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


kinematics_SP_data.distance.SP  = norm(kinematics_SP_data.position);
kinematics_SP_data              = distance_PL1___Sr___Sr(kinematics_SP_data);


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%