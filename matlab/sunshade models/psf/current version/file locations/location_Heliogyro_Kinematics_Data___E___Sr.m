
function    [excel_file] = ...
			........................................................................................................................................................ 
            location_Heliogyro_Kinematics_Data___E___Sr 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Define the location of the Excel file containing heliogyro positions and normal vectors.								                                %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	INPUT: Define the filename, folder, and complete location.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


excel_file_name                    = 'constellation_B2_11726crafts_2026-04-02.mat';
excel_file_folder                  = location_Excel_Folder___0___St;
excel_file.location                = fullfile(excel_file_folder, excel_file_name);
excel_file.name                    = excel_file_name;
excel_file.format                  = 'mat_v2';        % 'mat_v2' for V2 .mat files | 'xlsx_v1' for V1 .xlsx files


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Define the ranges within the Excel file that contain the position vectors, normal vectors, and initial time.
%	(Used only when format = 'xlsx_v1'. Retained here for backward compatibility with V1 constellations.)
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


excel_file.ranges.time_initial     = 'D8';
excel_file.ranges.position_vectors = 'A11:C20000';
excel_file.ranges.normal_vectors   = 'D11:F20000';


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%